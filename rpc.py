#!/usr/bin/env python2
import sys
import json
import numpy as np
def check_symmetric(mat):
	values = {}
	for row, indices in enumerate(mat['index']):
		for i, col in enumerate(indices):
			mat[row, col] = mat['value'][row][i]		
	for (row, col), val in values.iteritems():
		if mat[col, row] != val:
			print "not symmetric", row, col, val, mat[col, row]
def max_abs(*args):
	print "max_abs",
	for arg in args:
		if isinstance(arg, basestring):
			print arg,
		else:
			print np.absolute(arg).max(),
	print
def print_sum(name, mat):
	print "sum", name, np.sum(mat)
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
