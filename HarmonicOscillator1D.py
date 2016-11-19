# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum
# bond distance for n H2 atoms

import unittest
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.optimize import fmin

# Constants from LJ potential
# e: H well depth
# d: H van der Waals radius
_e = 7.607e-19    # kJ/mol
_d = 1.2          # angstroms

def generateAtoms(num):
	# Generates a random 3D array of n x 3 size
	array = np.random.rand(num, 3)
	array_size = array.shape
	return array, array_size


def euclideanDist(point1, point2):
	return spatial.distance.euclidean(point1, point2)


def calculateV(array):
	# Calculation of LJ potential using possible radii
	# Returns potential in kJ/mol as an list with closest neighbor(s)
	points = tuple(map(tuple, array))
	radius = []
	V = []

	# Query tree for nearest neighbors using KD tree
	# Returns indices of NN
	for n in points:
		output = [i for i in points if i != n]
		tree = spatial.KDTree(output)
		index = tree.query_ball_point(x=n, r=3)
		for i in index:
			radius.append(euclideanDist(tree.data[i], n))
	
	# Calculate potential for NNs
	for r in radius:
		V.append((4 * _e) * (((_d / r) ** 12) - ((_d / r) ** 6)))
	return V

	plt.plot(points)
	plt.show()


def moveMolecule(array):
	x0 = [0, 0, 0]
	minV = fmin(calculateV,x0)
	print minV

	for x,y,z in array:
		addArray = np.random.randint(low=-2, high=2, size=(len(array),3))
		while calculateV(array) != minV:
			array = np.add(array, addArray)


def putInPymol(array):
	# First cast array as strings
	def makeStr(myArray):
		return map(str, myArray)
	str_arr = map(makeStr, array)

	# Get dimensions of array and flatten to make life easier
	x = np.array(str_arr)
	flatArray = x.flatten()
	size = flatArray.shape

	# Put into PyMol readable format
	# Must iterate every 3 as we've flattened the array
	# Want format: ATOM      0  H0   U 0  10      x y z
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

	data2 = [([2, 3],
			 [1, 4],
			 [3, 2]),
			([2, 3],
			 [2, 5],
			 [4, 2],
			 [5, 1],
			 [4, 3])]

	easy = [([0, 0, 0],
			 [1, 1, 0],
			 [3, 1, 0])]


	# Test to see if generateAtoms generates
	# the right size array
	def test_generateAtoms(self):
		for [test_num, expected_size] in self.data1:
		    (actual_array, actual_size) = generateAtoms(test_num)
		    self.assertEqual(actual_size, expected_size)


	# Test to see if calculateV generates values
	# Compare to calculated values
	def test_calculateV(self):
		for array in self.easy:
			calculated_V = calculateV(array)
			actual_V = [-7.1181416371322853e-19, 
			-7.1181416371322853e-19, -1.3534136350801918e-19, 
			-1.3534136350801918e-19]
			self.assertEqual(calculated_V, actual_V)



	# Test putting arrays into Pymol coords
	def test_putInPymol(self):
   		for array in self.easy:
   			putInPymol(array)

   	
   	def test_moveMolecules(self):
   		for array in self.easy:
   			moveMolecule(array)
   	




if __name__ == "__main__":
	unittest.main()