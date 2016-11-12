# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum
# bond distance for n H2 atoms

import unittest
import cmath
import numpy as np
from scipy import spatial
import heapq

# constants from LJ potential
# e: H well depth
# d: H van der Waals radius
_e = 7.607e-19    # kJ/mol
_d = 1.2          # angstroms


def generateAtoms(num):
	# generates a random 3D array of n x 3 size
	array = np.random.rand(num, 3)
	array_size = array.shape
	return array, array_size

def calculateV(array):
	# calculation of LJ potential using possible radii
	# returns potential in kJ/mol
	npArray = np.array(array)

	# first generate matrix of nearest neighbors using KD tree
	# generate grid
	x, y, z = np.mgrid[0:5, 0:5, 0:5]

	# zip into tuples
	data = zip(x.ravel(), y.ravel(), z.ravel())

	# index nearest neighbors
	KDtree = spatial.KDTree(data)

	# then find all points within distance r of point x
	ball = KDtree.query_ball_point(npArray, 3)
	distance, index = KDtree.query(npArray)
	pts = npArray
	KDtree.query(pts)

	def euclideanDistance(x,y):
		dist = np.sqrt(sum([(a-b) ** 2 for (a,b) in zip(x,y)]))
		return dist

	# Utilize heap q algorithm (binary tree)
	# Sorts dataset to find smallest radius
	closestPoints = heapq.nsmallest(1, enumerate(pts), 
		key=lambda y: euclideanDistance(x, y[1]).any())

	print closestPoints
	radius = 

	# HALP
	_V = ((4 * _e) * (((_d / radius) ** 12) - ((_d / radius) ** 6)))
	return _V

	# used this site for help:
	# http://code.activestate.com/recipes/578434-a-simple-kd-tree-example-with-custom-euclidean-dis/

	

	"""
	# First I tried it my way :(
	x = np.array(array)
	flatArray = x.flatten()
	size = flatArray.shape
	print flatArray
	
	for n in range(size[0] - 5):
		radius_xy = ((float(flatArray[n]) - float(flatArray[n + 3])) ** 2 +
					 (float(flatArray[n + 1]) - float(flatArray[n + 4])) ** 2 +
					 (float(flatArray[n + 2]) - float(flatArray[n + 5])) ** 2)
	for n in range(size[0] - 5):
		radius_yz =  ((float(flatArray[n]) - float(flatArray[n + 3])) ** 2 +
					 (float(flatArray[n + 1]) - float(flatArray[n + 4])) ** 2 +
					 (float(flatArray[n + 2]) - float(flatArray[n + 5])) ** 2)

	for n in range(size[0] - 5):
		radius_xz =  ((float(flatArray[n]) - float(flatArray[n + 3])) ** 2 +
					 (float(flatArray[n + 1]) - float(flatArray[n + 4])) ** 2 +
					 (float(flatArray[n + 2]) - float(flatArray[n + 5])) ** 2)

		radius = cmath.sqrt(radius_xy + radius_yz + radius_xz)
		radiusList = []
		radiusList.append(radius)
		print radiusLi
	"""

def putInPymol(array):
	# First cast array as strings
	def makeStr(myArray):
		return map(str, myArray)
	str_arr = map(makeStr, array)

	# get dimensions of array and flatten to make life easier
	x = np.array(str_arr)
	flatArray = x.flatten()
	size = flatArray.shape

	# put into PyMol readable format
	# must iterate every 3 as we've flattened the array
	# want format: ATOM      0  H0   U 0  10      x y z
	# LATER FIX TO PUT M LOOP ON OUTSIDE--DEFINE ATOMS CORRECTLY
	file = open('testfile.pdb', 'w')
	for n in range(0, size[0] - 2, 3):
		file.write ('ATOM      ' + str(n) + '  H' + str(n) 
					+ '   U 0  10      ' )
		file.write(flatArray[n] + '.000  ' 
				 + flatArray[n + 1] + '.000  '
				 + flatArray[n + 2] + '.000  ')
		file.write('\n')
	file.close()



class Test(unittest.TestCase):
	data1 = [(2, (2, 3)), (3, (3, 3)), (0, (0, 3))]

	data2 = [([2, 3, 4],
			 [1, 4, 1],
			 [3, 2, 1]),
			([0, 0, 0],
			 [1, 1, 1],
			 [2, 2, 2]),
			([2, 3, 3],
			 [2, 5, 3],
			 [4, 2, 2],
			 [5, 1, 5],
			 [4, 2, 3])]

	data3 = [([0, 0, 0],
			 [1, 1, 1],
			 [2, 2, 2])]

	# test to see if generateAtoms generates
	# the right size array
	def test_generateAtoms(self):
		for [test_num, expected_size] in self.data1:
		    (actual_array, actual_size) = generateAtoms(test_num)
		    self.assertEqual(actual_size, expected_size)

	# test to see if calculateV generates values
	# COMPARE TO CALCULATED?
	def test_calculateV(self):
		for array in self.data3:
			calculateV(array)

	# test putting arrays into Pymol coords
	def test_putInPymol(self):
   		for array in self.data3:
   			putInPymol(array)
	



if __name__ == "__main__":
	unittest.main()