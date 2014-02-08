#!/usr/bin/env python2
import sys
import json
def check_symmetric(mat):
	values = {}
	for row, indices in enumerate(mat['index']):
		for i, col in enumerate(indices):
			mat[row, col] = mat['value'][row][i]		
	for (row, col), val in mat.iteritems():
		if mat[col, row] != val:
			print "not symmetric", row, col, val, mat[col, row]
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
	globals()['method'](*params)
