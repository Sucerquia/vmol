import vmol.vpython_with_img as vp


class VisualAngle:
    """New vpython object created out of triangles (semicircle). This allows to
    change the magnitude, orientation and position of the angle even after
    created.
    """
    def __init__(self, a: list, b: list, n: int,
                 origin: list = None, canvas: vp.canvas = None,
                 color: list = None, factor: float = 1):
        """
        Parameters
        ==========
        a: list
            vector side of the angle.
        b: list
            vector side of the angle.
        n: int
            number of intermedia triangles to define the arc.
        origin: list
            position of the vertex of the angle.
        canvas: vpython.canvas
            vpython scene to be added.
        color: list
            RGB color of the arc.
        factor: float
            factor of the shortest side to define the radius of the arc.
        """
        # number of intermedia vectors
        if color is None:
            self.color = vp.vector(0.5, 0.5, 0.5)
        else:
            self.color = self._asvector(color)
        self.n = n + 1
        self.factor = factor
        self.vertexes = [vp.vertex(pos=vp.vector(0, 0, 0),
                                   color=self.color,
                                   canvas=canvas) for _ in range(self.n)]

        a = self._asvector(a)
        b = self._asvector(b)
        if origin is None:
            origin = vp.vector(0, 0, 0)

        self.origin = vp.vertex(pos=self._asvector(origin),
                                color=self.color,
                                canvas=canvas)
        self.update_vertexes(self.origin.pos, a, b)

        self.triangles = []
        for i in range(self.n - 1):
            self.triangles.append(vp.triangle(v0=self.origin,
                                              v1=self.vertexes[i],
                                              v2=self.vertexes[i + 1],
                                              color=self.color,
                                              canvas=canvas))

    def _asvector(self, arraylike):
        """Changes an array to a vpython vector

        Returns
        =======
        (vpython.vector) transformed array-like into vpython.vector.
        """
        if not isinstance(arraylike, vp.vector):
            arraylike = vp.vector(*arraylike)
        return arraylike

    def update_vertexes(self, origin: list = None, a: list = None,
                        b: list = None, color: list = None) -> list:
        """Up-to-dates color, origin or the position of the intermedia
        arc-points based in new side vectors.

        Parameters
        ==========
        origin:
            position of the vertex of the arc.
        a: list
            vector side of the angle.
        b: list
            vector side of the angle.
        color: list
            RGB color of the arc.

        Returns
        =======
        (list) list of updated vertexes.
        """
        if a is None:
            a = self.vertexes[0].pos
        else:
            a = self._asvector(a)

        if b is None:
            b = self.vertexes[-1].pos
        else:
            b = self._asvector(b)

        if origin is not None:
            self.origin.pos = self._asvector(origin)
        if color is not None:
            color = self._asvector(color)
            self.origin.color = color
            for vertex in self.vertexes:
                vertex.color = color

        lena = a.mag
        lenb = b.mag
        axis = vp.cross(a, b).hat  # normal vector
        if axis.mag == 0:
            axis = vp.cross(a, vp.vector(0, 1, 0))
        lend = self.factor * min(lena, lenb)  # minumim length

        theta_total = vp.diff_angle(a, b)

        for i in range(self.n):
            theta = i * theta_total / (self.n - 1)
            self.vertexes[i].pos = self.origin.pos + \
                vp.rotate(lend*a.hat, angle=theta, axis=axis)
        return self.vertexes

    def transport(self, v0: list) -> list:
        """translate the angle as a rigid body
        Parameters
        ==========
        v0: list
            vector that to translate the arc.

        Returns
        =======
        (list) list of updated vertexes.
        """
        v0 = self._asvector(v0)
        self.origin.pos += v0
        for vertex in self.vertexes:
            vertex.pos += v0
        return self.vertexes

    def rotate(self, angle, axis):
        # TODO : to be implemented
        return True

    def hide(self):
        """hide the angle

        Returns
        =======
        (VisualAngle) hiden angle.
        """
        for tri in self.triangles:
            tri.visible = False
        return self

    def show(self):
        """shows the angle in the scene.

        Returns
        =======
        (VisualAngle) showed anfle.
        """
        for tri in self.triangles:
            tri.visible = True
        return self
