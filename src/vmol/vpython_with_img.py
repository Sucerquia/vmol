"""
from vpython._notebook_helpers import _isnotebook
from vpython._vector_import_helper import (vector, mag, norm, cross, dot, adjust_up,
                                    adjust_axis, object_rotate)
from vpython.vpython import baseObj, list_to_vec, distant_light, Mouse, color, Camera, box, print_to_string, event_return, _wait
from vpython.rate_control import rate

from math import sqrt, tan, pi
from inspect import signature
"""
# Cythonize the encode machinery?
import colorsys
from vpython.rate_control import rate
import platform

from math import sqrt, tan, pi

import time

# vpython provides clock in its namespace
clock = time.perf_counter

import sys
from vpython import __version__, __gs_version__
from vpython._notebook_helpers import _isnotebook
from vpython._vector_import_helper import (vector, mag, norm, cross, dot, adjust_up,
                                    adjust_axis, object_rotate)
from vpython import *
from vpython.vpython import list_to_vec, print_to_string, _wait
from inspect import signature


class canvas(baseObj):
    selected = None
    hasmouse = None
    maxVertices = 4.2e9  ## 2^32

    def __init__(self, **args):
        baseObj._canvas_constructing = True
        if _isnotebook:
            if 'img' in list(args.keys()):
                img = args['img']
            else:
                img = ''
            
            if 'height' in list(args.keys()):
                height = args['height']
            else:
                # because of this I just removed the
                # automatic creation of a canvas in __init__
                # just because it would mean a white unnecessary space 
                height = 500
                args['height'] = height
            
            # percentaje to of the 3d_canvas in the visualization
            if 'portion' in list(args.keys()):
                portion = args['portion']
            else:
                portion = 80

            from IPython.display import display, HTML, Javascript
            
            if img != '':
                display(HTML('<body>'
                            f'<div style="width: 100%; height: {height}px;">'
                            f'<div style="width: {portion}%;" id="glowscript" class="glowscript"></div>'
                            f'<div style="float: right; width: {100-portion}%;"> '
                            f'<img src="{img}"'
                            'style="object-fit:fill;'
                            f'height:{height}px;"/>'
                            ' </div>'
                            '</div>'
                            '</body>'))
                display(Javascript("""if (typeof Jupyter !== "undefined") { window.__context = { glowscript_container: $("#glowscript").removeAttr("id")};}else{ element.textContent = ' ';}"""))
            else:    
                display(HTML("""<div id="glowscript" class="glowscript"></div>"""))
                display(Javascript("""if (typeof Jupyter !== "undefined") { window.__context = { glowscript_container: $("#glowscript").removeAttr("id")};}else{ element.textContent = ' ';}"""))

        super(canvas, self).__init__()   ## get idx, attrsupdt

        self._constructing = True
        canvas.selected = self

        if 'lights' in args:
            raise AttributeError("Lights for a canvas can be assigned only after the canvas has been created.")

        self._objz = set()
        self.vertexCount = 0
        self._visible = True
        self._background = vector(0,0,0)
        self._ambient = vector(0.2, 0.2, 0.2)
        self._height = 400 # to match the GlowScript default
        self._width = 640
        self._align = 'none'
        self._fov = pi/3
        self._resizable = True

        # The following determine the view:
        self._range = 1 # user can alter with zoom
        self._axis = vector(0,0,-1) # user can alter with spin
        self._forward = vector(0,0,-1) # self.axis is primal internally; self._forward is now a synonym
        self._up = vector(0,1,0) # user with touch screen can rotate around z
        self._autoscale = True # set False if user zooms
        self._center = vector(0,0,0) # cannot be altered by user
        # Reject JavaScript canvas_update user values immediately following Python setting of values:
        self._set_range = False
        self._set_forward = False
        self._set_center = False
        self._set_up = False
        self._set_autoscale = False

        self._userzoom = True
        self._userspin = True
        self._userpan = True
        self._pixel_to_world = 0
        self._title = ''
        self._caption = ''
        self._mouse = Mouse(self)
        self._binds = {'mousedown':[], 'mouseup':[], 'mousemove':[],'click':[],
                        'mouseenter':[], 'mouseleave':[], 'keydown':[], 'keyup':[],
                        'redraw':[], 'draw_complete':[], 'resize':[]}
                        #'_compound':[]}
            # no key events unless notebook command mode can be disabled
        self._camera = Camera(self)
        self.title_anchor   = [self.idx, 1]  ## used by buttons etc.
        self.caption_anchor = [self.idx, 2]
        cmd = {"cmd": "canvas", "idx": self.idx}

    # send only nondefault values to GlowScript

        canvasVecAttrs = ['background', 'ambient','forward','up', 'center']
        canvasNonVecAttrs = ['visible', 'height', 'width', 'title','fov', 'range','align',
                             'autoscale', 'userzoom', 'userspin', 'userpan', 'title', 'caption']

        for a in canvasNonVecAttrs:
            if a in args:
                if args[a] is not None:
                    setattr(self, '_'+a, args[a])
                    cmd[a]= args[a]
                del args[a]

        for a in canvasVecAttrs:
            if a in args:
                aval = args[a]
                if not isinstance(aval, vector):
                    raise TypeError(a, 'must be a vector')
                setattr(self, '_'+a, vector(aval))
                cmd[a] = aval.value
                del args[a]

    # set values of user-defined attributes
        for key, value in args.items(): # Assign all other properties
            setattr(self, key, value)

        self._forward.on_change = self._on_forward_change
        self._up.on_change = self._on_up_change
        self._center.on_change = self._on_center_change

        self.appendcmd(cmd)
        self._constructing = False

        self._camera.follow = self.follow

        self.lights = [] # delete all lights created by glowcomm.js
        # Add the standard lighting (these lights will be added to self._lights):
        distant_light(direction=vector( 0.22,  0.44,  0.88), color=color.gray(0.8),
                      canvas=self)
        distant_light(direction=vector(-0.88, -0.22, -0.44), color=color.gray(0.3),
                      canvas=self)
        baseObj._canvas_constructing = False

    def follow(self, obj):
        if obj is None:
            self.addmethod('follow', 'None')
        elif callable(obj):
            b = box(visible=False)
            baseObj.follow_objects.append([b, obj, vector(1.2e15,3.4e14,-5.6e13)])
            self.addmethod('follow', b.idx)
        else:
            self.addmethod('follow', obj.idx)

    def select(self):
        canvas.selected = self
        self.addmethod('select','None')

    @classmethod
    def get_selected(cls):
        return cls.selected

    def delete(self):
        self.addmethod('delete','None')

    @property
    def title(self):
        return self._title
    @title.setter
    def title(self,value):
        self._title = value
        if not self._constructing:
            self.appendcmd({'title':value})

    @property
    def caption(self):
        return self._caption
    @caption.setter
    def caption(self,value):
        self._caption = value
        if not self._constructing:
            self.appendcmd({'caption':value})

    def append_to_title(self, *args):
        t = print_to_string(*args)
        self._title += t
        if not self._constructing:
            self.appendcmd({'append_to_title':t})

    def append_to_caption(self, *args):
        t = print_to_string(*args)
        self._caption += t
        if not self._constructing:
            self.appendcmd({'append_to_caption':t})

    @property
    def mouse(self):
        return self._mouse
    @mouse.setter
    def mouse(self,value):
        raise AttributeError('Cannot set scene.mouse')

    @property
    def camera(self):
        return self._camera
    @camera.setter
    def camera(self,value):
        raise AttributeError('Cannot set scene.camera')

    @property
    def visible(self):
        return self._visible
    @visible.setter
    def visible(self,value):
        self._visible = value
        if not self._constructing:
            self.addattr('visible')

    @property
    def resizable(self):
        return self._resizable
    @resizable.setter
    def resizable(self,value):
        self._resizable = value
        if not self._constructing:
            self.addattr('resizable')

    @property
    def background(self):
        return self._background
    @background.setter
    def background(self,value):
        self._background = value
        if not self._constructing:
            self.appendcmd({"background":value.value})

    @property
    def ambient(self):
        return self._ambient
    @ambient.setter
    def ambient(self,value):
        self._ambient = vector(value)
        if not self._constructing:
            self.addattr('ambient')

    @property
    def width(self):
        return self._width
    @width.setter
    def width(self,value):
        self._width = value
        if not self._constructing:
            self.appendcmd({"width":value})

    @property
    def height(self):
        return self._height
    @height.setter
    def height(self,value):
        self._height = value
        if not self._constructing:
            self.appendcmd({"height":value})

    @property
    def align(self): return self._align
    @align.setter
    def align(self,val):
        if not (val == 'left' or val == 'right' or val == 'none'):
            raise NameError("align must be 'left', 'right', or 'none' (the default).")
        self._align = val
        self.appendcmd({"align":val})

    @property
    def center(self):
        return self._center
    @center.setter
    def center(self,value):
        if isinstance(value, vector):
            self._center = self._set_center = vector(value)
            if not self._constructing:
                self.appendcmd({"center":value.value})
        else:
            raise TypeError('center must be a vector')

    @property
    def axis(self):
        return self._axis
    @axis.setter
    def axis(self,value):
        self._axis = self._set_forward = vector(value)
        if not self._constructing:
            self.appendcmd({"forward":value.value})

    @property
    def forward(self): # scene.forward is an external synonym for scene.axis
        return self._axis
    @forward.setter
    def forward(self,value):
        self._axis = self._set_forward = vector(value)
        if not self._constructing:
            self.appendcmd({"forward":value.value})

    @property
    def range(self):
        return self._range
    @range.setter
    def range(self,value):
        self._range = self._set_range = value
        if not self._constructing:
            self.appendcmd({"range":value})

    @property
    def up(self):
        return self._up
    @up.setter
    def up(self,value):
        self._up = self._set_up = value
        if not self._constructing:
            self.appendcmd({"up":value.value})

    @property
    def autoscale(self):
        return self._autoscale
    @autoscale.setter
    def autoscale(self,value):
        self._autoscale = self._set_autoscale = value
        if not self._constructing:
            self.appendcmd({"autoscale":value})

    @property
    def fov(self):
        return self._fov
    @fov.setter
    def fov(self,value):
        self._fov = value
        if not self._constructing:
            self.appendcmd({"fov":value})

    @property
    def userzoom(self):
        return self._userzoom
    @userzoom.setter
    def userzoom(self,value):
        self._userzoom = value
        if not self._constructing:
            self.appendcmd({"userzoom":value})

    @property
    def userspin(self):
        return self._userspin
    @userspin.setter
    def userspin(self,value):
        self._userspin = value
        if not self._constructing:
            self.appendcmd({"userspin":value})

    @property
    def userpan(self):
        return self._userpan
    @userpan.setter
    def userpan(self,value):
        self._userpan = value
        if not self._constructing:
            self.appendcmd({"userpan":value})

    @property
    def lights(self):
        return self._lights
    @lights.setter
    def lights(self, value):
        if type(value) is list and len(value) == 0:
            # JSON doesn't like an empty list
            self._lights = []
            self.appendcmd({"lights":'empty_list'}) # don't encode this unusual statement
        else:
            raise AttributeError("canvas.lights can be set only to [].")

    @property
    def pixel_to_world(self):
        return self._pixel_to_world
    @pixel_to_world.setter
    def pixel_to_world(self, value):
        raise AttributeError('pixel_to_world is read-only')

    def capture(self, filename, capture_labels=True):
        if not isinstance(filename, str):
            raise TypeError("'filename' for Capture must be a string.")

        if filename.endswith(".png"):
            filename += ".png"

        include_labels = "T" if capture_labels else "F"
        self.addmethod("capture", include_labels + filename)

    @property
    def objects(self):
        obs = []
        for ob in self._objz:
            if ob.visible:
                obs.append(ob)
        return obs
    @objects.setter
    def objects(self, *args1, **args ):
        raise AttributeError('objects is read-only')

    def objz(self, obj, operation):
        try:
            ii = (obj.idx > 0)  ## invalid object will not have .idx attribute
            if operation == 'add':
                self._objz.add(obj)
            elif operation == 'delete':
                self._objz.remove(obj)
        except:
            raise TypeError(obj + ' is not an object belonging to a canvas')

## key events conflict with notebook command mode; not permitted for now
    def handle_event(self, evt):  ## events and scene info updates
        global keysdownlist
        ev = evt['event']
        if ev == 'pick':
            self.mouse.setpick( evt )
            self._waitfor = True # what pick is looking for
        elif ev == '_compound': # compound, text, extrusion
            obj = self._compound
            p = evt['pos']
            if obj._objName == 'text':
                obj._length = p[0]
                obj._descender = p[1]
                obj._up.value = list_to_vec(evt['up'])
            else:
                if obj._objName == 'extrusion':
                    obj._color = obj._firstcolor # use the first segment color to represent multicolor extrusions
                # Set attribute_vector.value, which avoids nullifying the
                # on_change functions that detect changes in e.g. obj.pos.y
                obj._pos.value = list_to_vec(p)
                obj._origin = obj._pos
                s = evt['size']
                obj._size.value = obj._size0 = list_to_vec(s)
                obj._axis.value = obj._size._x*norm(obj._axis)
                obj._up.value = list_to_vec(evt['up'])
            self._waitfor = True # what compound and text and extrusion are looking for in _wait()
        elif ev == 'resize':
            if self.resizable and ('resize' in self._binds):
                self.width = evt['width']
                self.height = evt['height']
                del evt['width']
                del evt['height']
                for fct in self._binds['resize']:
                    # inspect the bound function and see what it's expecting
                    a = signature(fct)
                    if str(a) != '()':
                        fct(evt)
                    else:
                        fct()

        else: # pause/waitfor, update_canvas
            if 'pos' in evt:
                pos = evt['pos']
                evt['pos'] = list_to_vec(pos)
                self.mouse._pos = evt['pos']
            if 'ray' in evt:
                ray = evt['ray']
                evt['ray'] = list_to_vec(ray)
                self.mouse._ray = evt['ray']
            canvas.hasmouse = self
            if ev != 'update_canvas':   ## mouse events bound to functions, and pause/waitfor
                evt['canvas'] = self
                if ev[:3] != 'key':  # not a key event
                    self.mouse._alt = evt['alt']
                    self.mouse._shift = evt['shift']
                    self.mouse._ctrl = evt['ctrl']
                evt1 = event_return(evt)  ## turn it into an object
                for fct in self._binds[ev]:
                    # inspect the bound function and see what it's expecting
                    a = signature(fct)
                    if str(a) != '()':
                        fct( evt1 )
                    else:
                        fct()

                self._waitfor = evt1 # what pause and waitfor are looking for
            else:  ## user can change forward (spin), range/autoscale (zoom), up (touch), center (pan)
                if 'forward' in evt and self.userspin and not self._set_forward:
                    fwd = evt['forward']
                    self._axis = list_to_vec(fwd) # the fundamental meaning of scene.forward is scene.axis
                self._set_forward = False
                if 'up' in evt and self.userspin and not self._set_up:
                    cup = evt['up']
                    self._up = list_to_vec(cup)
                self._set_up = False
                if 'center' in evt and self.userpan and not self._set_center:
                    center = evt['center']
                    self._center = list_to_vec(center)
                self._set_center = False
                if 'range' in evt and self.userzoom and not self._set_range:
                    self._range = evt['range']
                self._set_range = False
                if 'autoscale' in evt and self.userzoom and not self._set_autoscale:
                    self._autoscale = evt['autoscale']
                self._set_autoscale = False
                if 'keysdown' in evt: keysdownlist = evt['keysdown']


    def bind(self, eventtype, whattodo):
        evts = eventtype.split()
        for evt in evts:
            if evt in self._binds:
                self._binds[evt].append(whattodo)
            else:
                raise TypeError(evt + ' is an illegal event type')
        self.addmethod('bind', eventtype)

    def unbind(self, eventtype, whatnottodo):
        evts = eventtype.split()
        for evt in evts:
            if evt in self._binds and whatnottodo in self._binds[evt]:
                self._binds[evt].remove(whatnottodo)
        self.addmethod('unbind', eventtype)

    def waitfor(self, eventtype):
        if 'textures' in eventtype: # textures are local; little need to wait
            eventtype = eventtype.replace('textures', '')
            if eventtype == '': return
        evts = ['redraw', 'draw_complete'] # wait for a render
        if eventtype in evts:
            baseObj.sent = False
            while baseObj.sent is False:
                rate(60)
        else:
            self.addmethod('waitfor', eventtype)
            _wait(self)
            return self._waitfor

    def pause(self,*s):
        if len(s) > 0:
            s = s[0]
            self.addmethod('pause', s)
        else:
            self.addmethod('pause', '')
        _wait(self)
        return self._waitfor

    def _on_forward_change(self):
        self.addattr('forward')

    def _on_up_change(self):
        self.addattr('up')

    def _on_center_change(self):
        self.addattr('center')

    def _ipython_display_(self): # don't print something when making an (anonymous) canvas
        pass
