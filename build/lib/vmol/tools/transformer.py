import numpy as np

# Alignment
class AtomicTrans:
    def rot_x(self, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[1, 0, 0],
                      [0, c, -s],
                      [0, s, c]])
        return R

    def rot_y(self, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, 0, s],
                      [0, 1, 0],
                      [-s, 0, c]])
        return R

    def rot_z(self, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector):
        """Apply the necessary rotations to  align a
        vector with the positive x axis
        """
        xyproj = vector.copy()
        xyproj[2] = 0
        phi = np.arcsin(vector[2] / np.linalg.norm(vector))
        theta = np.arccos(vector[0] / np.linalg.norm(xyproj))
        if vector[1] < 0:
            theta *= -1
        trans = np.dot(self.rot_y(phi), self.rot_z(-theta))
        return trans

    def align_plane(self, vector):
        """
        Rotation around x axis to set a vector in the xy plane
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1] / np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, atoms, trans, indexes=None):
        """Apply a transformation to all vector positions of the
        atoms object
        """
        if indexes is None:
            indexes = list(range(len(atoms)))

        new_positions = []
        for i, atom in enumerate(atoms):
            if i in indexes:
                new_positions.append(np.dot(trans, atom.position))
            else:
                new_positions.append(atom.position)
        atoms.set_positions(new_positions)

        return new_positions

    def xy_alignment(self, atoms, index1, index2, index3):
        """Transforme the positions of the atoms such that
        the atoms of indexes 1 and 2 are aligned in the
        x axis

        the atom 3 is in the xy plane"""
        # center
        center = (atoms[index1].position + atoms[index2].position) / 2
        atoms.set_positions(atoms.positions - center)
        axis = atoms[index2].position
        self.apply_trans(atoms, self.align_axis(axis))
        third = atoms[index3].position
        self.apply_trans(atoms, self.align_plane(third))

        return atoms
