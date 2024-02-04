from ase import Atoms
from ase.data.colors import jmol_colors
from vmol.tools.transformer import AtomicTrans
from vmol.tools.dofs import VisualAngle
from vpython import vector
import vpython as vp
import numpy as np

# removes the canvas by default and allows to create a new
# scene when a new MolView is created
vp.scene.delete()

class VMolecule(AtomicTrans):
    def __init__(self, atoms: Atoms,
                 color_scheme: dict = jmol_colors,
                 show_axis: bool = False,
                 alignment: list = None,
                 frame=0,
                 radius=0.2,
                 **kwargs) -> None:
        if 'background' in kwargs.keys():
            kwargs['background'] = self._asvector(kwargs['background'])

        # next line imports the creates
        self.scene = vp.canvas(**kwargs)

        # atoms and trajectory
        self.frame = frame
        if isinstance(atoms, Atoms):
            if alignment is not None:
                index1, index2, index3 = alignment
                atoms = self.xy_alignment(atoms, index1, index2, index3)
            self.trajectory = [atoms.copy()]
            self.frame = 0
        else:
            assert isinstance(atoms[0], Atoms), "atoms must be an ase.Atoms object or a trajectory of ase.Atoms"
            if alignment is not None:
                index1, index2, index3 = alignment
                atoms = [self.xy_alignment(config, index1, index2, index3)
                         for config in atoms]
            self.trajectory = [config.copy() for config in atoms]

        self.atoms = atoms[self.frame]
        self.axis = None
        self.show_axis(show=show_axis)
        self.color_scheme = color_scheme
        self.vatoms = []
        self.dofs = {}
        # Unfortunately, vpython has problems removing objects (check vpython.baseObj.delete,
        # then I just keep them invisible, in case to reload,
        # just change to visibles dofs
        self.removed_dofs = {}
        # TODO: add option to define atoms radius
        self.add_atoms(radius=radius)
        # vpython such that the user can use vpython commands as well
        self.master = vp
        
    # region 1_vpython-tools
    def _asvector(self, arraylike):
        if not isinstance(arraylike, vp.vector):
            arraylike = vp.vector(*arraylike)
        return arraylike
    
    
    def update_obj(self, obj, **kwargs):
        for attr, val in kwargs.items():
            # The next conditional assumes that none of the vpyhon
            # objects has lists as parameters.
            if isinstance(val, (np.ndarray,
                                list,
                                tuple)):
                val = vector(*val)
            setattr(obj, attr, val)
        return obj
    
    def setting_canvas(self, **kwargs):
        scene = self.update_obj(self.scene, **kwargs)
        return scene
    
    def show_axis(self, show: bool = False,
                  size: float = 1.0):
        if self.axis is None:
            x = vp.arrow(pos=vector(0,0,0),
                         axis=vector(size,0,0),
                         color=vector(1,0,0),
                         canvas=self.scene)
            y = vp.arrow(pos=vector(0,0,0),
                         axis=vector(0,size,0),
                         color=vector(0,1,0),
                         canvas=self.scene)
            z = vp.arrow(pos=vector(0,0,0),
                         axis=vector(0,0,size),
                         color=vector(0,0,1),
                         canvas=self.scene)
            self.axis = [x, y, z]
        for ax in self.axis:
            ax.visible = show
        
        return self.axis

    def download_image(self, name):
        return self.scene.capture(name)
    # endregion
    
    # region 2_Atoms
    def add_atoms(self, atoms: Atoms = None,
                  color_scheme: np.ndarray=jmol_colors,
                  radius: float = 0.2):
        if atoms is None:
            atoms = self.atoms
        for i, atom in enumerate(atoms):
            sphere = vp.sphere(pos=vector(*atom.position),
                               color=vector(*color_scheme[atom.number]),
                               canvas=self.scene,
                               radius=radius)
            sphere.index = i + 1
            self.vatoms.append(sphere)
    
    def update_atom(self, index, **kwargs):
        for attr, val in kwargs.items():
            # The next conditional assumes that none of the vpyhon
            # objects has lists as parameters.
            if isinstance(val, (np.ndarray,
                                list,
                                tuple)):
                val = vector(*val)
            setattr(self.vatoms[index], attr, val)
        return self.vatoms[index]
    # endregion
    
    # region 2_DOFs_1_bonds
    def add_bond(self, atom1index, atom2index,
                 color=None, radius=0.2):
        """Add a bond between two atoms:
        atom1 and atom2

        Parameters
        ==========
        atom1index (and atom2index): int
            Indexes of the atoms to be connected according with g09
            convention.

        color: list. Default gray([0.5, 0.5, 0.5])
            RGB triplet.

        radius: float. Default 0.1
            Radius of the bond.

        Output
        ======

        Return the bonds in the system
        
        """
        bond_extreme = self.vatoms[atom1index - 1].pos
        bond_axis = self.vatoms[atom2index - 1].pos - \
                    self.vatoms[atom1index - 1].pos
        
        if color is None:
            # add ball+stick colors as nglview
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)

        indexes = [atom1index, atom2index]
        if atom1index > atom2index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)
        
        # If already exists
        if name in self.dofs.keys():
            self.update_obj(self.dofs[name],
                            pos=bond_extreme,
                            axis=bond_axis,
                            radius=radius,
                            color=color)
            return self.dofs[name]
        # If it was previously removed. It means it is there but it is invisible.
        elif name in self.removed_dofs.keys():
            self.dofs[name] = self.removed_dofs[name]
            del self.removed_dofs[name]
            self.update_obj(self.dofs[name],
                            pos=bond_extreme,
                            axis=bond_axis,
                            visible=True,
                            radius=radius,
                            color=color)
            return self.dofs[name]

        # In case it was not created previously:
        b = vp.cylinder(pos=bond_extreme,
                        axis=bond_axis,
                        color=color,
                        radius=radius)
        
        b.indices = np.array([atom1index, atom2index, 0, 0])

        self.dofs[name] = b

        return self.dofs[name]
    
    def add_bonds(self, atoms1indexes, atoms2indexes, colors=None, radii: float=0.07):
        """Add a bond between each pair of atoms corresponding to
        two lists of atoms:
        atoms1 and atoms.

        Parameters
        ==========
        atom1index (and atom2index): int
            Indexes of the atoms to be connected according with g09
            convention.
        colors: list of color lists. Default all gray([0.5, 0.5, 0.5])
            RGB triplets for each of the bonds. It can be one a triplet
            in case of just one color in all bonds.
        radii: float or list of floats. Default 0.1
            radius of each bond.

        Output
        ======

        Return the bonds in the system
        """

        if colors is None:
            colors = [0.5, 0.5, 0.5]

        if not isinstance(colors[0], (list, np.ndarray)):
            colors = [colors] * len(atoms1indexes)

        if type(radii) is not list:
            radii = [radii] * len(atoms1indexes)

        assert len(atoms1indexes) == len(atoms2indexes), \
            "The number of atoms in both lists must be the same"
        assert len(atoms1indexes) == len(colors), \
            "The number of colors in must be the same as the number of atoms"
        assert len(atoms1indexes) == len(radii), \
            "The number of radii must be the same as the number of atoms"

        for i in range(len(atoms1indexes)):
            self.add_bond(atoms1indexes[i],
                          atoms2indexes[i],
                          colors[i],
                          radii[i])
        return self.bonds
    
    def remove_bond(self, atom1index, atom2index):
        """Remove a bond between two atoms:
        atoms1 and atoms2.

        Parameters
        ==========

        atom1index (and atom2index): int
            Indexes of the atoms that are connected. This bond
            will be removed.

        Output
        ======

        Return the bonds in the system
        """
        indexes = [atom1index, atom2index]
        indexes.sort()
        name = ''.join(str(i).zfill(3) for i in indexes)
        if name in self.dofs.keys():
            self.dofs[name].visible = False
            self.removed_dofs[name] = self.dofs[name]
            del self.dofs[name]
        return self.dofs
    
    def remove_bonds(self, atoms1indexes=None, atoms2indexes=None):
        """Remove several bonds in the plot between two list of atoms:
        atoms1 and atoms2.

        Parameters
        ==========

        atom1index (and atom2index): list[int]
            Indexes of the atoms that are connected.

        Note: if atoms2 is None, all bonds with atoms1 will me removed.
        If atoms1 and atoms2 are None, all bonds in the structure are
        removed.
        """

        if (atoms1indexes is None) and (atoms2indexes is None):
            for name in self.dofs.keys():
                if len(name) == 6:
                    self.dofs[name].visible = False
                    self.removed_dofs[name] = self.dofs[name]
                    del self.dofs[name]
            return self.removed_dofs

        elif (atoms1indexes is not None) and (atoms2indexes is None):
            for name in self.dofs.keys():
                for index in atoms1indexes:
                    if (str(index).zfill(3) in name) and (len(name) == 6):
                        self.dofs[name].visible = False
                        self.removed_dofs[name] = self.dofs[name]
                        del self.dofs[name]
            return self.removed_dofs

        else:
            assert len(atoms1indexes) == len(atoms2indexes), \
                "The number of atoms in both lists must be the same"
            for index1, index2 in zip(atoms1indexes, atoms2indexes):
                self.remove_bond(index1, index2)
            return self.removed_dofs

    def remove_all_bonds(self):
        """ remove all bonds"""
        return self.remove_bonds()
    # endregion
    
    # region 2_DOFs_2_angles
    def add_angle(self, atom1index, atom2index, atom3index,
                  color=None, n=20, factor=0.7):
        """Add an angle to between three atoms:
        atom1, atom2 and atom3
        - with the vertex in the atom2

        Parameters
        ==========

        atom1index, atom2index and atom3index: int
            Indexes of the three atoms that defines the angle.
        color: color list. Default all gray([0.5, 0.5, 0.5])
            RGB triplet.
        n: int. Default 10
            number of intermedia points to add in the arc of
            the angle.

        Output
        ======
        Return the angles in the system
        """
        if color is None:
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)

        #if self.is_trajectory:
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
        # If it was previously removed. It means it is there but it is invisible.
        elif name in self.removed_dofs.keys():
            self.dofs[name] = self.removed_dofs[name]
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            self.dofs[name].show()
            del self.removed_dofs[name]
            return self.dofs[name]

        # In case it was not created previously:
        self.dofs[name] = VisualAngle(side1, side2, n,
                                      origin=vertex, color=color,
                                      scene=self.scene,
                                      factor=factor)

        self.dofs[name].indices = np.array([atom1index, atom2index,
                                            atom3index, 0])
        
        return self.dofs[name]

    def remove_angle(self, atom1index, atom2index, atom3index):
        indexes = [atom1index, atom2index, atom3index]
        if atom1index > atom3index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.dofs.keys():
            self.dofs[name].hide()
            self.removed_dofs[name] = self.dofs[name]
            del self.dofs[name]
        return self.removed_dofs
    
    def remove_all_angles(self):
        for dof in self.dofs.keys():
            if len(dof) == 9:
                self.dofs[dof].hide()
                self.removed_dofs[dof] = self.dofs[dof]
                del self.dofs[dof]
        return self.removed_dofs
    # endregion
    
    # region 2_DOFs_3_dihedrals
    def add_dihedral(self, atom1index, atom2index, atom3index,
                     atom4index, color=None, n=20, factor=0.7):
        """Add an dihedral angle between four atoms:
        atom1, atom2, atom3 and atom4
        - with the vertex in the midle of the atom 2 and 3

        Parameters
        ==========

        atom1index, atom2index, atom3index and atom4index: int
            Indexes of the three atoms that defines the angle.
        color: color list. Default all gray([0.5, 0.5, 0.5])
            RGB triplet.
        n: int. Default 10
            number of intermedia points to add in the arc of
            the angle.

        Output
        ======
        Return the dihedral angles
        """
        if color is None:
            color = vp.vector(0.5, 0.5, 0.5)
        color = self._asvector(color)

        #if self.is_trajectory:
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
        # If it was previously removed. It means it is still there but it is invisible.
        elif name in self.removed_dofs.keys():
            self.dofs[name] = self.removed_dofs[name]
            self.dofs[name].update_vertexes(origin=vertex,
                                            a=side1,
                                            b=side2,
                                            color=color)
            self.dofs[name].show()
            del self.removed_dofs[name]
            return self.dofs[name]

        # In case it was not created previously:
        self.dofs[name] = VisualAngle(side1, side2, n,
                                      origin=vertex,
                                      color=color,
                                      scene=self.scene,
                                      factor=factor)
        self.dofs[name].indices = np.array([atom1index, atom2index,
                                            atom3index, atom4index])
        return self.dofs[name]
    
    def remove_dihedral(self, atom1index, atom2index, atom3index, atom4index):
        indexes = [atom1index, atom2index, atom3index, atom4index]
        if atom1index > atom4index:
            indexes = indexes[::-1]
        name = ''.join(str(i).zfill(3) for i in indexes)

        if name in self.dofs.keys():
            self.dofs[name].hide()
            self.removed_dofs[name] = self.dofs[name]
            del self.dofs[name]
        return self.removed_dofs

    def remove_all_dihedrals(self):
        for dof in self.dofs.keys():
            if len(dof) == 12:
                self.dofs[dof].hide()
                self.removed_dofs[dof] = self.dofs[dof]
                del self.dofs[dof]
        return self.removed_dofs
    # endregion
    
    # region 2_DOFs_4_arbitrary
    def add_dof(self, dof, **kwargs):
        """Add the degree of freedom to the molecule image

        Parameters
        ==========

        dof: tuple
            label of the degree of freedom according with g09 convention.
        
        kwargs:
            arguments corresponding with the DOF you want to add.

        Example
        =======
            i=(1, 2) means a bond between atoms 1 and 2
            i=(1, 2, 3) means an angle between atoms 1, 2 and 3
            i=(1, 2, 3, 4) means a dihedral angle between atoms 1, 2, 3 and 4
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

    # region update_frame
    def update_frame(self, frame):
        self.frame = frame
        self.atoms = self.trajectory[self.frame]

        for i, atom in enumerate(self.vatoms):
            atom.pos = vp.vector(*self.atoms[i].position)
        for dof in self.dofs.values():
            self.add_dof(dof.indices)
        return frame
    # endregion

    
