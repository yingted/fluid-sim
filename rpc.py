#!/usr/bin/env python2
import sys
import json
import numpy as np
import thread
import threading
import traceback
import collections
import OpenGL
import fcntl
import termios
import struct
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
def mat_to_array(mat):
	index = mat['index']
	value = mat['value']
	ret = np.zeros((len(index), max(map(max, index))+1))
	for row, indices in enumerate(index):
		for i, col in enumerate(indices):
			ret[row, col] = value[row][i]
	return ret
def check_symmetric(mat):
	for row, col in zip(*np.nonzero(mat != mat.transpose())):
		if row < col:
			print "not symmetric", row, col, mat[row, col], mat[col, row]
def max_abs(*args):
	print "max_abs",
	for arg in args:
		if isinstance(arg, basestring):
			print arg,
		else:
			print np.absolute(arg).max(),
	print
def deciles(name, mat):
	print "deciles", name,
	for decile in xrange(0, 101, 10):
		print "% 5.04f"%np.percentile(mat, decile),
	print
def print_sum(name, mat):
	print "sum", name, np.sum(mat)
def opengl(w, h, frameskip=False):
	def decorate(redraw):
		lock = threading.Lock()
		q = collections.deque(maxlen = 1 if frameskip else None)
		dirty = [True]
		def idle():
			with lock:
				if not dirty[0]:
					q.popleft()
					dirty[0] = True
				if q:
					glutPostRedisplay()
				elif not frameskip:
					sys.exit()
		def display():
			with lock:
				ret._frame, args = q[0]
				dirty[0] = False
			try:
				redraw(*args)
			except:
				traceback.print_exc()
				thread.interrupt_main()
				thread.exit()
		def ret(*args):
			lock.acquire()
			dirty[0] = True
			if ret.window is not None:
				q.append((q[-1][0]+1, args))
				lock.release()
				return
			q.append((0, args))
			def glut_thread():
				glutInit()
				glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
				glutInitWindowSize(w, h)
				ret.window = glutCreateWindow(redraw.__name__)
				glutIdleFunc(idle)
				glutDisplayFunc(display)
				lock.release()
				glutMainLoop()
			thread = threading.Thread(target=glut_thread)
			thread.daemon = frameskip
			thread.start()
		ret._state = -1, None
		ret.window = None
		ret.__name__ = redraw.__name__
		return ret
	return decorate
@opengl(500, 500, frameskip=False)
def draw(dx, dy, phi, mx, my, bx, by, r):
	phi = np.array(phi)
	h, w = phi.shape
	phi = ((np.concatenate((
		phi.transpose(),
		[[0]*((-w)%4)]*h,
	), axis=1)/r+1)*128).clip(0, 255).astype("uint8").tostring()
	glPushMatrix()
	if True:

		glColorMask(1, 1, 1, 1)
		glClearColor(0., 0., 0., 0.)
		glClear(GL_COLOR_BUFFER_BIT)

		glColor3f(1, 1, 1)
		glColorMask(1, 1, 1, 1)

		#glShadeModel(GL_SMOOTH)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST) # instead of GL_LINEAR
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, w, h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, phi)
		glEnable(GL_TEXTURE_2D)
		glBegin(GL_QUADS)
		if True:
			glTexCoord2f(0., 0.)
			glVertex2f(-1., -1.)
			glTexCoord2f(1., 0.)
			glVertex2f( 1., -1.)
			glTexCoord2f(1., 1.)
			glVertex2f( 1.,  1.)
			glTexCoord2f(0., 1.)
			glVertex2f(-1.,  1.)
		glEnd()
		glDisable(GL_TEXTURE_2D)

		glTranslatef(-1., -1., 0.)
		glScalef(2./w, 2./h, 1.)

		theta = np.linspace(0, 2*np.pi, num=np.ceil(2*np.pi/np.arccos(1-1/(500./w*r))))
		circle = (np.array([np.cos(theta), np.sin(theta)])*r).transpose()
		glVertexPointerf(circle)
		glEnableClientState(GL_VERTEX_ARRAY)

		m = np.array([mx, my]).transpose()
		m = np.concatenate((m[:1], np.diff(m, axis=0)))
		for r, color in ((500./w+1)/(500./w), (0, 0, 1, 1)), (1., (1, 0, 0, 1)):
			glColorMask(*color)
			glColor3f(*color)
			glPushMatrix()
			glScalef(r, r, r)
			for x, y in m/r:
				glTranslatef(x, y, 0)
				glDrawArrays(GL_TRIANGLE_FAN, 0, len(circle))
			glPopMatrix()

		glDisableClientState(GL_VERTEX_ARRAY)

		glColor3f(0, 0, 0)
		glColorMask(1, 1, 1, 1)
		boundary = np.array([bx, by]).transpose()
		glVertexPointerf(boundary)
		glEnableClientState(GL_VERTEX_ARRAY)
		glDrawArrays(GL_POINTS, 0, len(boundary))
		glDisableClientState(GL_VERTEX_ARRAY)

		glColor3f(0, 1, 0)
		dx = np.array(dx)
		dy = np.array(dy)
		pos = np.concatenate((
			np.mgrid[map(slice, dx.shape)].reshape(2, dx.size)+np.transpose([[0, .5]]*dx.size),
			np.mgrid[map(slice, dy.shape)].reshape(2, dy.size)+np.transpose([[.5, 0]]*dy.size),
		), axis=1)
		vel = np.concatenate((
			np.array([dx, np.zeros(dx.shape)]).reshape(2, dx.size),
			np.array([np.zeros(dy.shape), dy]).reshape(2, dy.size),
		), axis=1)
		vectors = np.rollaxis(np.array((pos, pos+vel)), -1)
		glVertexPointerf(vectors)
		glEnableClientState(GL_VERTEX_ARRAY)
		glDrawArrays(GL_LINES, 0, 2*len(vectors))
		glDisableClientState(GL_VERTEX_ARRAY)

		glutSwapBuffers()
	glPopMatrix()
def convert_param(param):
	if isinstance(param, dict):
		if isinstance(param.get('index'), list) and isinstance(param.get('value'), list):
			return mat_to_array(param)
	return param
def echo(*args):
	for arg in args:
		print arg,
	print
def terminal_size():
	try:
		h, w, hp, wp = struct.unpack('HHHH', fcntl.ioctl(sys.stdout.fileno(), termios.TIOCGWINSZ, struct.pack('HHHH', 0, 0, 0, 0)))
	except IOError:
		h = w = None
	return h, w
def main():
	np.set_printoptions(precision=3, suppress=True, linewidth=terminal_size()[1])
	for line in sys.stdin:
		try:
			o = json.loads(line)
		except ValueError:
			continue
		try:
			method = o['method']
			params = o['params']
		except KeyError:
			continue
		globals()[method](*map(convert_param, params))
if __name__ == "__main__":
	main()
