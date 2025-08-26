from ase import Atoms
from ase.data.colors import jmol_colors
from vmol.tools.transformer import AtomicTrans
from vmol.tools.dofs import VisualAngle
from vpython import vector
import vpython as vp
import numpy as np
from IPython.display import clear_output
from ase.geometry.analysis import Analysis

# hides the canvas by default and allows to create a new
# scene when a new MolView is created
vp.scene.delete()


class VMolecule(AtomicTrans):
    """Tool to visualize molecules, equivalent to nglview or vmd"""
    def __init__(self, atoms: Atoms,
                 color_scheme: dict = jmol_colors,
                 show_axis: bool = False,
                 alignment: list = None,
                 frame=0,
                 radius=0.2,
                 default_bonds=False,
                 **kwargs) -> None:
        """
        Parameters
        ==========
        atoms: ase.Atoms
            Atoms object with the molecule structure
        color_scheme: dict
            color scheme of the atoms. Default = jml_colors
        show_axis: bool
            True to show xyz axis. Default = False
        alignment: list[int]
            list of three indexes corresponding to the indexes of the atoms in
            the xy plane. 1-based indexing. The first two atoms are set on the
            x axis.
        frame: int
            In case to be a trajectory, the index of the captured frame. In
            case of being a single point configuration, frame is always 0.
        radius: float
            set the size of the atoms. Default = 0.2
        **kwargs of Vpython canvas
        """
        self._hook = False
        if 'background' in kwargs.keys():
            kwargs['background'] = self._asvector(kwargs['background'])

        # vpython such that the user can use vpython commands as well
        self.vp = vp

        # next line imports the creates
        clear_output()
        self.scene = vp.canvas(**kwargs)
        self.create_caption()

        # atoms and trajectory
        self.frame = frame
        if isinstance(atoms, Atoms):
            if alignment is not None:
                index1, index2, index3 = alignment
                atoms = self.xy_alignment(atoms, index1, index2, index3)
            self.trajectory = [atoms]
            self.frame = 0
        else:
            assert isinstance(atoms[0], Atoms), "atoms must be an ase.Atoms" +\
                " object or a list of ase.Atoms"
            if alignment is not None:
                index1, index2, index3 = alignment
                atoms = [self.xy_alignment(config, index1, index2, index3)
                         for config in atoms]
            self.trajectory = [config.copy() for config in atoms]

        self.atoms = self.trajectory[self.frame]
        self.axis = None
        self.selected = None
        self.show_axis(show=show_axis)
        self.color_scheme = color_scheme
        self.vatoms = []
        self.dofs = {}
        # Unfortunately, vpython has problems removing objects
        # (check vpython.baseObj.delete,
        # then I just keep them invisible, in case to reload,
        # just change to visibles dofs
        self.hidden_objs = {}
        self.add_atoms(radius=radius)
        
        if default_bonds:
            self.add_def_bonds(self.atoms)

        self._hook = True

    def __getattribute__(self, name):
        attr = super().__getattribute__(name)
        update_frame = super().__getattribute__('update_frame')
        if callable(attr) and \
            name != 'update_frame' and \
            name != '_hook' and \
            super().__getattribute__('_hook'):
            def wrapper(*args, **kwargs):
                outp = attr(*args, **kwargs)
                update_frame()
                return outp  # Original method
            return wrapper
        else:
            return attr

    # region 1_vpython-tools
    def _asvector(self, arraylike) -> vector:
        """Changes an array to a vpython vector

        Returns
        =======
        (vpython.vector) transformed array-like into vpython.vector.
        """
        if not isinstance(arraylike, vp.vector):
            arraylike = vp.vector(*arraylike)
        return arraylike

    def update_obj(self, obj, **kwargs):
        """Changes attributes of an object. Specially useful to update
        positions, size and colors in vpython scene.

        Parameters
        ==========
        obj: obj
            object to be modified.
        **kwargs: attributes and values to be modified.

        Returns
        =======
        (obj) Modified Object

        Usage
        =====
        VMolecule.update_obj(vpython.sphere, radius=10)
        """
        for attr, val in kwargs.items():
            # The next conditional assumes that none of the vpyhon
            # objects has lists as parameters.
            if isinstance(val, (np.ndarray,
                                list,
                                tuple)):
                val = vector(*val)
            setattr(obj, attr, val)
        return obj

    def setting_canvas(self, **kwargs) -> vp.canvas:
        """Setting vpython scene.

        \*\*kwargs: vpython.canvas attibutes
        """
        scene = self.update_obj(self.scene, **kwargs)
        return scene
    
    def create_caption(self):

        # Note: when this function fails, does not return any message. Just
        # wtext stop working. cchanges texts and creating new ones.
        self.caption = vp.wtext(text='', canvas=self.scene)
        def identify():
            at = self.scene.mouse.pick
            if at is not None and hasattr(at, 'index'):
                self.selected = at
                self.caption.text = str(self.vatoms[at.index - 1].info)
                self.caption.text = self.caption.text.replace(',', '\n')
                self.caption.text = self.caption.text.replace('{', ' ')
                self.caption.text = self.caption.text.replace('}', '')
                self.caption.text = self.caption.text.replace("'", '')
            else:
                self.selected = None
                self.caption.text = ''
        self.scene.bind('click', identify)

    def show_axis(self, show: bool = True,
                  size: float = 1.0,
                  center: vector = vector(0, 0, 0)) -> list:
        """show xyz axis

        Parameters
        ==========
        show: bool
            True to show it. Default = True
        size: float
            thickness of the arrows
        center: vector or array
            position of the origin

        Returns
        =======
        (list) list with the xyz vpytho arrows
        """
        if self.axis is None:
            x = vp.arrow(pos=center,
                         axis=vector(size, 0, 0),
                         color=vector(1, 0, 0),
                         canvas=self.scene)
            y = vp.arrow(pos=center,
                         axis=vector(0, size, 0),
                         color=vector(0, 1, 0),
                         canvas=self.scene)
            z = vp.arrow(pos=center,
                         axis=vector(0, 0, size),
                         color=vector(0, 0, 1),
                         canvas=self.scene)
            self.axis = [x, y, z]
        for ax in self.axis:
            ax.visible = show

        return self.axis

    def download_image(self, name):
        """Save displayed scene in a png file.

        Parameters
        ==========
        name: str
            name of the image file.

        Note: So far, I couldn't avoid to render the saving-file interface that
        promps, even when a name is given. This makes annoyng to render
        n-images for an animation."""
        return self.scene.capture(name)
    # endregion

        
    def add_image_to_canvas(self, fig):
        """
        Add the a matplotlib figure to the canvas in a notebook (other
        environments are not working properly yet).

        Parameters
        ==========
        fig: plt.Figure
            figure to be added to the canvas.

        Return
        ======
        (io.BytesIO) buffer that contains the image.

        Note
        ----
        This function assumes that the height of the figure is the same
        as the height of the canvas. Always take into account the that the size
        of the canvas is in pixels and the size of the matplotlib figure is in
        inches.

        Warning
        -------
        For some strange reason, this function does not work properly if a
        vpython canvas hasn't been rendered before.
        """
        if self.vp._notebook_helpers._isnotebook:
            import io
            import base64
            from IPython.display import display, HTML, Javascript

            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=fig.get_dpi())
            buf.seek(0)

            # Convert the PNG to a base64 string
            img_base64 = base64.b64encode(buf.read()).decode('utf-8')
            buf.close()

            # custom container with a "slot" for GlowScript + the image
            display(HTML(f"""
            <div id="glowscript" style="display: flex">
                <div id="vpython_target" style="flex: 0;"></div>
                <div>
                    <img src="data:image/png;base64,{img_base64}" 
                         style="height:{self.scene._height}px;"/>
                </div>
            </div>
            """))

            # tell GlowScript/VPython to render into the left slot
            display(Javascript("""
            if (typeof Jupyter !== "undefined") {
                window.__context = { glowscript_container:
                               $("#vpython_target").removeAttr("id") };
            } else {
                element.textContent = ' ';
            }
            """))

            return buf
        else:
            raise EnvironmentError("add_image_to_canvas only works in a " + 
                                   "Jupyter notebook environment.")

    # region 2_Atoms
    def add_atoms(self, atoms: Atoms = None,
                  color_scheme: dict = jmol_colors,
                  radius: float = 0.2) -> list:
        """Add vpython spheres that are representations of atoms.

        Parameters
        ==========
        atoms: ase.Atoms
            set of atoms to be rendered.
        color_scheme: dict
            color scheme of the atoms. Default = jml_colors
        radius: float
            set the size of the atoms. Default = 0.2

        Returns
        =======
        (list) set of created vpython spheres, namely, atoms in the scene.
        """
        if atoms is None:
            atoms = self.atoms

        for i, atom in enumerate(atoms):
            sphere = vp.sphere(pos=vector(*atom.position),
                               color=vector(*color_scheme[atom.number]),
                               canvas=self.scene,
                               radius=radius)
            sphere.index = i + 1
            sphere.info = {'Index': i + 1,
                           'Atom': atom.symbol,
                           'Pos_x': atom.x,
                           'Pos_y': atom.y,
                           'Pos_z': atom.z,
                           'Charge': atom.charge}
            self.vatoms.append(sphere)

    def add_def_bonds(self, atoms, engine='ase'):
        """
        Add the default bounds predicted by ASE.

        Parameters
        ==========
        atoms: ase.Atoms
            atoms to be connected by default bonds.
        
        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        if engine == 'ase':
            ana = Analysis(atoms)
            unique_contact = ana.unique_bonds[0]
        elif engine == 'vmd':
            from ase.io import write
            from myutils.ase_utils.molecules import vmd_connectivity_unique
            import os


            write('tmp_Mview.xyz', atoms)
            unique_contact = vmd_connectivity_unique('tmp_Mview.xyz')
            os.remove('tmp_Mview.xyz')


        bonds1 = []
        bonds2 = []
        for i, contact in enumerate(unique_contact):
            for j in contact:
                bonds1.append(i + 1)
                bonds2.append(j + 1)
        self.add_bonds(bonds1, bonds2)

    def hide_atom(self, atomindex: int) -> vp.sphere:
        """Hide one atom in the scene.

        Parameters
        ==========
        atomindex: int
            Index of the atom to be hidden. 1-based indexing.

        Returns
        =======
        (vpython.sphere) hidden vpython atom
        """
        self.vatoms[atomindex - 1].visible = False
        return self.vatoms[atomindex - 1]

    def show_atom(self, atomindex: int) -> vp.sphere:
        """Show one atom in the scene.

        Parameters
        ==========
        atomindex: int
            Index of the atom to be shown. 1-based indexing.

        Returns
        =======
        (vpython.sphere) hidden vpython atom
        """
        self.vatoms[atomindex - 1].visible = True
        return self.vatoms[atomindex - 1]

    def update_atom(self, index: int, **kwargs) -> vp.sphere:
        """Update atom properties

        Parameters
        ==========
        index: int
            index of the atom to be modified. 1-based indexing.
        **kwargs: arguments to be modified

        Returns
        =======
        (vpython.sphere) sphere representing the atom
        """
        self._hook = False
        atom = self.update_obj(self.vatoms[index - 1], **kwargs)
        self.atoms[index - 1].position = [self.vatoms[index - 1].pos.x,
                                          self.vatoms[index - 1].pos.y,
                                          self.vatoms[index - 1].pos.z]
        self._hook = True

        return atom
    # endregion

    # region 2_DOFs_1_bonds
    def add_bond(self, atom1index: int, atom2index: int,
                 color: vector = None, radius: float = 0.1) -> vp.cylinder:
        """
        Add a bond between two atoms.

        Parameters
        ==========
        atom1index: int
            Index of atom to be connected. 1-based indexing.
        atom2index: int
            Index of atom to be connected. 1-based indexing.
        color: list. Default=None
            RGB triplet. If None, it will be gray.
        radius: float. Default 0.1
            Radius of the bond.

        Returns
        =======
        (vpython.cylinder) representation of the bond.
        """
        bond_extreme = self.vatoms[atom1index - 1].pos
        bond_axis = self.vatoms[atom2index - 1].pos - \
            self.vatoms[atom1index - 1].pos

        indexes = [atom1index, atom2index]
        if atom1index > atom2index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        # If already exists
        if name in self.dofs.keys():
            if color is None:
                color = self.dofs[name].color
            self.update_obj(self.dofs[name],
                            pos=bond_extreme,
                            axis=bond_axis,
                            radius=radius,
                            color=color)
            return self.dofs[name]
        # If it was previously hidden. It means it is there but it is
        # invisible.
        elif name in self.hidden_objs.keys():
            if color is None:
                color = self.hidden_objs[name].color    
            self.dofs[name] = self.hidden_objs[name]
            del self.hidden_objs[name]
            self.update_obj(self.dofs[name],
                            pos=bond_extreme,
                            axis=bond_axis,
                            visible=True,
                            radius=radius,
                            color=color)
            return self.dofs[name]
        
        if color is None:
            # Add half color like in nglvew ball-sticks
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)

        # In case it was not created previously:
        b = vp.cylinder(pos=bond_extreme,
                        axis=bond_axis,
                        color=color,
                        radius=radius,
                        canvas=self.scene)

        b.indices = np.array([atom1index, atom2index, 0, 0])

        self.dofs[name] = b

        return self.dofs[name]

    def add_bonds(self, atoms1indexes: list,
                  atoms2indexes: list,
                  colors: list = None,
                  radii: float = 0.07) -> dict:
        """
        Add a bond between each pair of atoms corresponding to the ordered
        elements two lists of atoms.

        Parameters
        ==========
        atoms1indexes: list[int]
            Indexes of the atoms to be connected. 1-based indexing.
        atoms1indexes: list[int]
            Indexes of the atoms to be connected. 1-based indexing.
        colors: color-lists. Default=[0.5, 0.5, 0.5]
            RGB triplets for each of the bonds. It can be one a triplet
            in case of just one color in all bonds.
        radius: float or list of floats. Default=0.1
            radius of each bond.

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """

        if colors is None:
            colors = [0.5, 0.5, 0.5]

        if not isinstance(colors[0], (list, np.ndarray)):
            colors = [colors] * len(atoms1indexes)

        if type(radii) is not list:
            radii = [radii] * len(atoms1indexes)

        assert len(atoms1indexes) == len(atoms2indexes), \
            "The two lists defining the bonds have to have the same " + \
            "number of indexes"
        assert len(atoms1indexes) == len(colors), \
            "The number of colors in must be the same as the number of atoms"
        assert len(atoms1indexes) == len(radii), \
            "The number of radii must be the same as the number of atoms"

        for i in range(len(atoms1indexes)):
            self.add_bond(atoms1indexes[i],
                          atoms2indexes[i],
                          colors[i],
                          radii[i])
        return self.dofs

    def hide_bond(self, atom1index: int, atom2index: int) -> dict:
        """Hide one bond between two atoms:
        atoms1 and atoms2.

        Parameters
        ==========
        atom1: int
            Index of one of the atoms that are connected. This bond will be
            hidden. 1-based indexing.
        atom2: int
            Index of one of the atoms that are connected. This bond will be
            hidden. 1-based indexing.

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        indexes = [atom1index, atom2index]
        indexes.sort()
        name = ''.join(str(i).zfill(3) for i in indexes)
        if name in self.dofs.keys():
            self.dofs[name].visible = False
            self.hidden_objs[name] = self.dofs[name]
            del self.dofs[name]
        return self.dofs

    def hide_bonds(self, atoms1indexes: int = None,
                     atoms2indexes: int = None) -> dict:
        """Hide several bonds in the plot between two list of atoms:
        atoms1 and atoms2.

        Parameters
        ==========
        atoms1indexes: list[int]
            Indexes of the atoms that are connected. 1-based indexing.
        atoms2indexes: list[int]
            Indexes of the atoms that are connected. 1-based indexing.

        Returns
        =======
        (dict) all hidden DOFs in the system. keys -> dof names,
        values -> vpython.objects

        Note: if atoms2 is None, all bonds with atoms1 will me hidden.
        If atoms1 and atoms2 are None, all bonds in the structure are
        hidden.
        """

        if (atoms1indexes is None) and (atoms2indexes is None):
            todel = []
            for name in self.dofs.keys():
                if len(name) == 6:
                    self.dofs[name].visible = False
                    self.hidden_objs[name] = self.dofs[name]
                    todel.append(name)
            for name in todel:
                del self.dofs[name]

            return self.hidden_objs

        elif (atoms1indexes is not None) and (atoms2indexes is None):
            for name in self.dofs.keys():
                for index in atoms1indexes:
                    if (str(index).zfill(3) in name) and (len(name) == 6):
                        self.dofs[name].visible = False
                        self.hidden_objs[name] = self.dofs[name]
                        del self.dofs[name]
            return self.hidden_objs

        else:
            assert len(atoms1indexes) == len(atoms2indexes), \
                "The number of atoms in both lists must be the same"
            for index1, index2 in zip(atoms1indexes, atoms2indexes):
                self.hide_bond(index1, index2)
            return self.hidden_objs

    def hide_all_bonds(self):
        """Hide all bonds

        Returns
        =======
        (dict) all hidden DOFs in the system. keys -> dof names,
        values -> vpython.objects"""
        return self.hide_bonds()
    # endregion

    # region 2_DOFs_2_angles
    def add_angle(self, atom1index: int,
                  atom2index: int,
                  atom3index: int,
                  color: list = None,
                  n: int = 20,
                  factor: float = 0.7) -> dict:
        """Add an angle to between three atoms with the vertex in the atom2.

        Parameters
        ==========
        atom1index: int
            Index of the first atom that defines the angle. 1-based indexing.
        atom2index: int
            Index of the second atom that defines the angle. 1-based indexing.
        atom3index: int
            Index of the third atom that defines the angle. 1-based indexing.
        color: color list. Default all gray([0.5, 0.5, 0.5])
            RGB triplet.
        n: int. Default 10
            number of intermedia points to add in the arc of
            the angle.
        factor:
            factor of the minimim bond lenght to define the size of the arc.

        Returns
        =======
        (vmol.tools.dofs.VisualAngle) Just created angle.
        """
        if color is None:
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)
        # TODO: hide next block
        # if self.is_trajectory:
        #    atoms = self.atoms[self.viewer.view.frame]

        indexes = [atom1index, atom2index, atom3index]
        if atom1index > atom3index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        vertex = self.vatoms[atom2index - 1].pos
        side1 = self.vatoms[atom1index - 1].pos - vertex
        side2 = self.vatoms[atom3index - 1].pos - vertex

        # If already exists.
        if name in self.dofs.keys():
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            return self.dofs[name]
        # If it was previously hidden. It means it is there but it is
        # invisible.
        elif name in self.hidden_objs.keys():
            self.dofs[name] = self.hidden_objs[name]
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            self.dofs[name].show()
            del self.hidden_objs[name]
            return self.dofs[name]

        # In case it was not created previously:
        self.dofs[name] = VisualAngle(side1, side2, n,
                                      origin=vertex, color=color,
                                      canvas=self.scene,
                                      factor=factor)

        self.dofs[name].indices = np.array([atom1index, atom2index,
                                            atom3index, 0])

        return self.dofs[name]

    def hide_angle(self, atom1index: int,
                     atom2index: int,
                     atom3index: int) -> dict:
        """Hide one angle between three atoms:
        atoms1, atoms2 and atom3.

        Parameters
        ==========
        atom1index: int
            Index of the first atom that defines the angle. 1-based indexing.
        atom2index: int
            Index of the second atom that defines the angle. 1-based indexing.
        atom3index: int
            Index of the third atom that defines the angle. 1-based indexing.

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        indexes = [atom1index, atom2index, atom3index]
        if atom1index > atom3index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.dofs.keys():
            self.dofs[name].hide()
            self.hidden_objs[name] = self.dofs[name]
            del self.dofs[name]
        return self.hidden_objs

    def hide_all_angles(self):
        """Hide all the angles in the scene

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        for dof in self.dofs.keys():
            if len(dof) == 9:
                self.dofs[dof].hide()
                self.hidden_objs[dof] = self.dofs[dof]
                del self.dofs[dof]
        return self.hidden_objs
    # endregion

    # region 2_DOFs_3_dihedrals
    def add_dihedral(self, atom1index: int,
                     atom2index: int, atom3index: int,
                     atom4index: int,
                     color: list = None,
                     n: int = 20,
                     factor: float = 0.7) -> VisualAngle:
        """Add a dihedral angle between four atoms with the vertex in the midle
        of the atom 2 and 3.

        Parameters
        ==========
        atom1index: int
            Index of the first atom that defines the angle. 1-based indexing.
        atom2index: int
            Index of the second atom that defines the angle. 1-based indexing.
        atom3index: int
            Index of the third atom that defines the angle. 1-based indexing.
        atom4index: int
            Index of the fourth atom that defines the angle. 1-based indexing.
        color: color list
            RGB triplet. Default all gray([0.5, 0.5, 0.5])
        n: int. Default=20
            number of intermedia points to add in the arc of
            the angle.
        factor: float
            factor of the minimim bond lenght to define the size of the arc.

        Returns
        =======
        (vmol.tools.dofs.VisualAngle) Just created dihedral angle.
        """
        if color is None:
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)

        # TODO: delete next block
        # if self.is_trajectory:
        #    atoms = self.atoms[self.viewer.view.frame]

        indexes = [atom1index, atom2index, atom3index, atom4index]
        if atom1index > atom4index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        axis = (self.vatoms[atom3index - 1].pos -
                self.vatoms[atom2index - 1].pos)
        vertex = 0.5 * (self.vatoms[atom3index - 1].pos +
                        self.vatoms[atom2index - 1].pos)
        axis1 = (self.vatoms[atom1index - 1].pos -
                 self.vatoms[atom2index - 1].pos)
        axis2 = (self.vatoms[atom4index - 1].pos -
                 self.vatoms[atom3index - 1].pos)

        side1 = axis1 - vp.proj(axis1, axis)
        side2 = axis2 - vp.proj(axis2, axis)

        # If already exists.
        if name in self.dofs.keys():
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            return self.dofs[name]
        # If it was previously hidden. It means it is still there but it is
        # invisible.
        elif name in self.hidden_objs.keys():
            self.dofs[name] = self.hidden_objs[name]
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            self.dofs[name].show()
            del self.hidden_objs[name]
            return self.dofs[name]

        # In case it was not created previously:
        self.dofs[name] = VisualAngle(side1, side2, n,
                                      origin=vertex,
                                      color=color,
                                      canvas=self.scene,
                                      factor=factor)
        self.dofs[name].indices = np.array([atom1index, atom2index,
                                            atom3index, atom4index])
        return self.dofs[name]

    def hide_dihedral(self, atom1index: int, atom2index: int,
                        atom3index: int, atom4index: int):
        """Hide one dihedral angle defined by three atoms:
        atoms1, atoms2 and atom3.

        Parameters
        ==========
        atom1index: int
            Index of the first atom that defines the angle. 1-based indexing.
        atom2index: int
            Index of the second atom that defines the angle. 1-based indexing.
        atom3index: int
            Index of the third atom that defines the angle. 1-based indexing.
        atom4index: int
            Index of the fourth atom that defines the angle. 1-based indexing.

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        indexes = [atom1index, atom2index, atom3index, atom4index]
        if atom1index > atom4index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.dofs.keys():
            self.dofs[name].hide()
            self.hidden_objs[name] = self.dofs[name]
            del self.dofs[name]
        return self.hidden_objs

    def hide_all_dihedrals(self):
        """Hide all the dihedrals in the scene

        Returns
        =======
        (dict) all the DOFs in the system. keys -> dof names,
        values -> vpython.objects
        """
        for dof in self.dofs.keys():
            if len(dof) == 12:
                self.dofs[dof].hide()
                self.hidden_objs[dof] = self.dofs[dof]
                del self.dofs[dof]
        return self.hidden_objs
    # endregion

    # region 2_DOFs_4_arbitrary
    def add_dof(self, dof: list, **kwargs):
        """Add the degree of freedom to the scene

        Parameters
        ==========
        dof: list
            label of the degree of freedom according with g09 convention.
        **kwargs: arguments corresponding with the DOF you want to add.

        Returns
        =======
        Object used to represent the DOF. For further detail, check add_<dof>,
        where <dof> could be bond, angle or dihedral.

        Usage
        =======
        VMolecule.add_dof(dof), where:
            dof=(1, 2) means a bond between atoms 1 and 2
            dof=(1, 2, 3) means an angle between atoms 1, 2 and 3
            dof=(1, 2, 3, 4) means a dihedral angle between atoms 1, 2, 3 and 4
        """

        types = ["bond", "angle", "dihedral"]
        type_dof = types[np.count_nonzero(dof) - 2]

        if type_dof == "bond":
            index1 = dof[0]
            index2 = dof[1]
            return self.add_bond(index1, index2, **kwargs)

        elif type_dof == "angle":
            index1 = dof[0]
            index2 = dof[1]
            index3 = dof[2]
            return self.add_angle(index1, index2, index3, **kwargs)

        elif type_dof == "dihedral":
            index1 = dof[0]
            index2 = dof[1]
            index3 = dof[2]
            index4 = dof[3]
            return self.add_dihedral(index1, index2, index3,
                                     index4, **kwargs)
        else:
            raise TypeError(f"{dof} is not an accepted degree of freedom.")
    # endregion

    def seen_atoms(self):
        """
        Return an ase.Atoms object with the atoms that are displayed.
        """
        subset = []
        for i, atom in enumerate(self.vatoms):
            if atom.visible:
                self.atoms[i].position = [self.vatoms[i].pos.x,
                                          self.vatoms[i].pos.y,
                                          self.vatoms[i].pos.z]
                subset.append(i)

        return self.atoms[subset]
    
    def transform(self, trans, *args, **kwargs):
        """
        Applies a transformation to all the atoms in the trajectory.

        Parameters
        ==========
        trans: str
           name of the transformation you want to apply. For more information,
           visit vmol.tools.transformer
        *args and **kwargs of the transformation function.

        Return
        ======
        (list or Atoms) set of images in the trajectory. when there is only one
        structure, the return is a list with only one object.
        """
        for atoms in self.trajectory:
            getattr(self, trans)(atoms, *args, **kwargs)
        self.update_frame(self.frame)

        return self.trajectory

    # region update_frame
    def update_frame(self, frame=None):
        self._hook = False
        if frame is not None:
            self.frame = frame
        self.atoms = self.trajectory[self.frame]

        for i, atom in enumerate(self.vatoms):
            atom.pos = vp.vector(*self.atoms[i].position)
            atom.info['Pos_x'] = atom.pos.x
            atom.info['Pos_y'] = atom.pos.y
            atom.info['Pos_z'] = atom.pos.z
                           
        for dof in self.dofs.values():
            self.add_dof(dof.indices)
        self._hook = True
        return frame
    # endregion
