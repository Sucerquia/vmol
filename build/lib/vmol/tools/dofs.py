import vpython as vp


class VisualAngle:
    """New vpython object created out of triangles (semicircle).
    This allows to change the magnitude, orientation and position of the angle even
    after created.
    """
    def __init__(self, a, b, n, origin=None, scene=None, color=None,
                factor: float = 1):
        # number of intermedia vectors
        if color is None:
            self.color = vp.vector(0.5, 0.5, 0.5)
        else:
            self.color = self._asvector(color)
        self.n = n + 1
        self.factor = factor
        self.vertexes = [vp.vertex(pos=vp.vector(0,0,0),
                                   color=self.color) for _ in range(self.n)]

        a = self._asvector(a)
        b = self._asvector(b)
        if origin is None:
            origin = vp.vector(0, 0, 0)

        self.origin =  vp.vertex(pos=self._asvector(origin),
                                 color=self.color)
        self.update_vertexes(self.origin.pos, a, b)

        self.triangles = []
        for i in range(self.n - 1):
            self.triangles.append(vp.triangle(v0=self.origin,
                                              v1=self.vertexes[i],
                                              v2=self.vertexes[i + 1],
                                              color=self.color,
                                             canvas=scene))
    def _asvector(self, arraylike):
        if not isinstance(arraylike, vp.vector):
            arraylike = vp.vector(*arraylike)
        return arraylike

    def update_vertexes(self, origin=None, a=None, b=None, color=None):
        if a is None:
            a = self.vertexes[0].pos
        if b is None:
            b = self.vertexes[-1].pos
        if origin is not None:
            self.origin.pos = self._asvector(origin)
        if color is not None:
            color = self._asvector(color)
            self.origin.color = color
            for vertex in self.vertexes:
                vertex.color = color

        lena = a.mag
        lenb = b.mag
        axis = vp.cross(a, b).hat # normal vector
        if axis.mag == 0:
            axis = vp.cross(a, vp.vector(0,1,0))
        lend = self.factor * min(lena, lenb) # minumim length

        theta_total = vp.diff_angle(a, b)

        for i in range(self.n):
            theta = i * theta_total / (self.n - 1)
            self.vertexes[i].pos = self.origin.pos + vp.rotate(lend*a.hat, angle=theta, axis=axis)
        return self.vertexes
    
    def transport(self, v0):
        v0 = self._asvector(v0)
        self.origin.pos += v0
        for vertex in self.vertexes:
            vertex.pos += v0

    def rotate(self, angle, axis):
        #TODO : to be implemented
        return True
    
    def hide(self):
        for tri in self.triangles:
            tri.visible = False
        return self
    
    def show(self):
        for tri in self.triangles:
            tri.visible = True
        return self
