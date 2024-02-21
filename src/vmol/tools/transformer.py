import numpy as np
from ase import Atoms


class AtomicTrans:
    def rot_x(self, angle: float) -> np.ndarray:
        """Rotation matrix around x axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the x axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[1, 0, 0],
                      [0, c, -s],
                      [0, s, c]])
        return R

    def rot_y(self, angle: float) -> np.ndarray:
        """Rotation matrix around y axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the y axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, 0, s],
                      [0, 1, 0],
                      [-s, 0, c]])
        return R

    def rot_z(self, angle: float) -> np.ndarray:
        """Rotation matrix around z axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the z axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector: np.ndarray) -> np.ndarray:
        """Apply the necessary rotations to  align a vector with the positive x
        axis

        Parameters
        ==========
        vector: numpy array
            vector to be aligned with the positive x axis.

        Return
        ======
        (numpy.array) transformation matrix.
        """
        xyproj = vector.copy()
        xyproj[2] = 0
        phi = np.arcsin(vector[2] / np.linalg.norm(vector))
        theta = np.arccos(vector[0] / np.linalg.norm(xyproj))
        if vector[1] < 0:
            theta *= -1
        trans = np.dot(self.rot_y(phi), self.rot_z(-theta))
        return trans

    def align_plane(self, vector: np.ndarray) -> np.ndarray:
        """Rotation around x axis to set a vector in the xy plane

        Parameters
        ==========
        vector: numpy array
            vector to be aligned with the positive xy plane.

        Returns
        =======
        (numpy.array) rotation matrix around x axis lo align the vector with
        the xy plane.
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1] / np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, atoms: Atoms,
                    trans: np.ndarray,
                    shift: np.ndarray = None,
                    indexes: list = None) -> np.ndarray:
        """Apply a transformation to all vector positions of the atoms object

        Parameters
        ==========
        atoms: ase.Atoms
            Atoms object to apply the transformation.
        trans: numpy.array
            matrix transformation to apply to the atom positions.
        shift: numpy.array
            translation transformation.
        indexes: list
            indexes of the atoms to be transformed.

        Returns
        =======
        (list) new set of positions of all atoms. Even those that were not
        transformed.
        """
        if shift is None:
            shift = [0, 0, 0]
        if indexes is None:
            indexes = list(range(len(atoms)))

        new_positions = []
        for i, atom in enumerate(atoms):
            if i in indexes:
                new_positions.append(np.dot(trans, atom.position) + shift)
            else:
                new_positions.append(atom.position + shift)
        atoms.set_positions(new_positions)

        return new_positions

    def xy_alignment(self, atoms: Atoms, index1: int,
                     index2: int, index3: int = None,
                     center: int = None) -> Atoms:
        """Transform the positions of the atoms such that the atoms of indexes1
        and index2 are aligned in the x axis. The atom with index3 would be in
        the xy plane in case to be given.

        Parameters
        ==========
        index1 and index2: int
            indexes of the atoms to be aligned with the x-axis. The positive
            direction of x would go from atom 1 to atom 2.
        index 3: int (optional)
            The atom with index 3 would be in the xy plane in case to be given.
        Center: int (optional)
            It must be index1 or index2, that means the atom with this index
            will be placed in the origin. In case center=None (default), the
            origin would be in the geometrical center between atoms with index1
            and index2.

        Return
        ======
        (numpy.array)[#natoms x 3float] new xyz positions of the N atoms. It
            changes the positions of the internal atoms object.
        """
        # center
        if center == index1 or center == index2:
            center = atoms[center].position
        else:
            pos1 = atoms[index1].position
            pos2 = atoms[index2].position
            center = (pos1 + pos2) / 2

        atoms.set_positions(atoms.positions - center)
        axis = atoms[index2].position - atoms[index1].position
        self.apply_trans(atoms, self.align_axis(axis))
        if index3 is not None:
            third = atoms[index3].position
            self.apply_trans(atoms, self.align_plane(third))

        return atoms
