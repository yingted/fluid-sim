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
import itertools
import operator
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
def mat_to_array(mat):
	index = mat['index']
	value = mat['value']
	ret = np.zeros((len(index), max(itertools.chain((0,),itertools.chain.from_iterable(index)))+1))
	for row, indices in enumerate(index):
		for i, col in enumerate(indices):
			ret[row, col] = value[row][i]
	return ret
def _square(mat):
	if mat.shape[1] < mat.shape[0]:
		mat = np.concatenate((mat, [[0]*(mat.shape[0]-mat.shape[1])]*len(mat)), axis=1)
	return mat
def check_symmetric(mat):
	mat = _square(mat)
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
		parent = [None]
		q = collections.deque(maxlen = 1 if frameskip else None)
		dirty = [True]
		def idle():
			with lock:
				if not dirty[0] and len(q) > 1:
					q.popleft()
					dirty[0] = True
				if len(q)-(not dirty[0]):
					glutPostRedisplay()
				elif not frameskip and not parent[0].is_alive():
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
			parent[0] = threading.current_thread()
			thread = threading.Thread(target=glut_thread)
			thread.daemon = frameskip
			thread.start()
		ret._state = -1, None
		ret.window = None
		ret.__name__ = redraw.__name__
		return ret
	return decorate
winsize = 500
def _image(arr):
	arr = np.array(arr)
	h, w = arr.shape
	return h, w, ((np.concatenate((
		arr.transpose(),
		[[0]*((-w)%4)]*h,
	), axis=1)+1)*128).clip(0, 255).astype("uint8").tostring()
@opengl(winsize, winsize, frameskip=False)
def draw(solid_phi, dx, dy, phi, bx, by):
	h, w, phi = img_phi = _image(phi)
	_, _, solid_phi = img_solid_phi = _image(solid_phi)
	glPushMatrix()
	if True:

		glColorMask(1, 1, 1, 1)
		glClearColor(0., 0., 0., 0.)
		glClear(GL_COLOR_BUFFER_BIT)

		for (img_h, img_w, img), color in (img_solid_phi, (0, 0, 1)), (img_phi, (1, 1, 0)):
			glPushMatrix()
			glScalef(img_w/float(w), img_h/float(h), 1.)
			glColor3f(*color)
			glColorMask(*(color+(1,)))

			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST) # instead of GL_LINEAR
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
			glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, img_w, img_h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img)
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
			glPopMatrix()

		glTranslatef(-1., -1., 0.)
		glScalef(2./w, 2./h, 1.)

		glColor3f(0, 1, 0)
		glColorMask(1, 1, 1, 1)
		dx = np.array(dx)
		dy = np.array(dy)
		pos = np.concatenate((
			np.mgrid[map(slice, dx.shape)].reshape(2, dx.size)+np.transpose([[0, .5]]*dx.size),
			np.mgrid[map(slice, dy.shape)].reshape(2, dy.size)+np.transpose([[.5, 0]]*dy.size),
		), axis=1)
		vel = np.concatenate((
			np.array([dx, np.zeros(dx.shape)]).reshape(2, dx.size)*(float(winsize)/w),
			np.array([np.zeros(dy.shape), dy]).reshape(2, dy.size)*(float(winsize)/h),
		), axis=1)
		vectors = np.rollaxis(np.array((pos, pos+vel)), -1)
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf(vectors)
		glDrawArrays(GL_LINES, 0, 2*len(vectors))
		glDisableClientState(GL_VERTEX_ARRAY)

		glColor3f(1, 0, 0)
		glColorMask(1, 1, 1, 1)
		boundary = np.array([bx, by]).transpose()
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointerf(boundary)
		glDrawArrays(GL_POINTS, 0, len(boundary))
		glDisableClientState(GL_VERTEX_ARRAY)

		glutSwapBuffers()
	glPopMatrix()
@opengl(winsize, winsize, frameskip=False)
def draw_quad(cells):
	vecs = {k : map(operator.itemgetter(k), cells) for k in cells[0]}
	N = sum(vecs['leaf'])
	vecs = {k : np.fromiter(itertools.compress(v, vecs['leaf']), np.float, N) for k, v in vecs.iteritems()}
	del vecs['leaf']
	_draw_quad_vecs(**vecs)
def _draw_quad_vecs(solid_phi, phi, x, y, r):
	N = len(r)

	glColorMask(1, 1, 1, 1)
	glClearColor(0., 0., 0., 0.)
	glClear(GL_COLOR_BUFFER_BIT)

	glEnableClientState(GL_COLOR_ARRAY)
	glEnableClientState(GL_VERTEX_ARRAY)
	colors = np.array([[0, 0, 1]]*N)
	glColorPointer(3, GL_FLOAT, 0, colors)
	lx = x-r
	ly = y-r
	hx = x+r
	hy = y+r
	cells = np.array([lx, ly, hx, ly, hx, hy, lx, hy]).transpose().reshape((4*N, 2))
	glVertexPointerf(cells)
	glDrawArrays(GL_QUADS, 0, N)
	glDisableClientState(GL_VERTEX_ARRAY)
	glDisableClientState(GL_COLOR_ARRAY)

	glutSwapBuffers()
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
