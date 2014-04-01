#!/usr/bin/env python2
import sys
import json
import numpy as np
import thread
import threading
import traceback
import collections
import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
def check_symmetric(mat):
	values = {}
	for row, indices in enumerate(mat['index']):
		for i, col in enumerate(indices):
			values[row, col] = mat['value'][row][i]
	for (row, col), val in values.iteritems():
		if values.get((col, row)) != val:
			print "not symmetric", row, col, val, values.get((col, row))
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
def update_phi(phi, mx, my, r=1.02):
	h = len(phi)
	w = len(phi[0])
	phi = ((np.concatenate((
		np.array(phi).transpose(),
		[[0]*((-w)%4)]*h,
	), axis=1)/r+1)*128).clip(0, 255).astype("uint8").tostring()
	glPushMatrix()
	if True:

		glColorMask(1, 1, 1, 1)
		glClearColor(0., 0., 0., 0.)
		glClear(GL_COLOR_BUFFER_BIT)

		glColor4f(1, 1, 1, 1)
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

		theta = np.linspace(0, 2*np.pi, num=500./w*r*2*np.pi/7)
		circle = (np.array([np.cos(theta), np.sin(theta)])*r).transpose()
		glVertexPointerf(circle)
		glEnableClientState(GL_VERTEX_ARRAY)

		for r, color in ((500./w+1)/(500./w), (0, 0, 1, 1)), (1., (1, 0, 0, 1)):
			glColorMask(*color)
			glColor3f(*color)
			for x, y in zip(mx, my):
				glPushMatrix()
				glTranslatef(x, y, 0)
				glScalef(r, r, r)
				glDrawArrays(GL_TRIANGLE_FAN, 0, len(circle))
				glPopMatrix()

		glDisableClientState(GL_VERTEX_ARRAY)

		glutSwapBuffers()
	glPopMatrix()
def main():
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
		globals()[method](*params)
if __name__ == "__main__":
	main()
