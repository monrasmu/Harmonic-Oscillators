# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum
# bond distance for n H2 atoms
# ONLY WORKING FOR 3 ATOMS RIGHT NOW

import unittest
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
import matplotlib.lines as line
from mpl_toolkits.mplot3d import Axes3D

# Constants from LJ potential
# e: H well depth
# d: H van der Waals radius
# Wilhelm, E.; Battino, R. J Chem Physics 1971, 55, 4012
_e = 0.09794  # kj/mol
_d = 1.09     # angstroms

# Number of atoms
# Variable to be changed as needed
num = 3

def doItAll(num):
	# It does it all!
	# Generates random atoms and moves them to minimums
	array = generateAtoms(num)
	moveMolecule(array)


def generateAtoms(num):
	# Generates a random 3D array of num x 3 size
	array = np.random.rand(num, 3)
	array_size = array.shape
	array = np.squeeze(array)
	return array


def euclideanDist(point1, point2):
	# Returns euclidean distance between two points
	return spatial.distance.euclidean(point1, point2)


def NNSearch(array):
	# Query tree for nearest neighbors using KD tree
	# Returns radius as list of NN distances
	points = tuple(map(tuple, array))
	radius = []

	for n in points:
		output = [i for i in points if i != n]
		tree = spatial.KDTree(output)
		index = tree.query_ball_point(x=n, r=3)
		for i in index:
			radius.append(euclideanDist(tree.data[i], n))
	return radius

def calculateV(array):
	# Calculation of LJ potential using possible radii
	# Returns potential in kJ/mol as an list with closest neighbor(s)
	V = []
	
	radius  = NNSearch(array)
	
	# Calculate potential for NNs
	for r in radius:
		V.append((4 * _e) * (((_d / r) ** 12) - ((_d / r) ** 6)))
	return V


def moveMolecule(array):
	# Moves molecules around to find minimum bond distance
	# Stop when potential below LJ potential for H2
	# Plots path of points
	# Puts final output into Pymol

	points = []
	x = []
	y = []
	z = []

	# Moves molecules in random directions by adding random array
	# Cutoff value chosen from minimum of function
	while all(V > -0.09 for V in calculateV(array)):
		addArray = np.random.random_sample(size=(len(array),3))
		print addArray
		array = np.add(array, addArray)
		points.append(array)

	# Extract points for plotting
	points = np.squeeze(points)

	for i in points:
		for j in i:
			x.append(j[0])
			y.append(j[1])
			z.append(j[2])

	finalPoint1 = (x[-3], y[-3], z[-3])
	finalPoint2 = (x[-2], y[-2], z[-2])
	finalPoint3 = (x[-1], y[-1], z[-1])
	finalPoints = []
	finalPoints.append(finalPoint1)
	finalPoints.append(finalPoint2)
	finalPoints.append(finalPoint3)

	putInPymol(finalPoints)

	print 'between 1 and 2: ', euclideanDist(finalPoint1, finalPoint2)
	print 'between 1 and 3: ', euclideanDist(finalPoint1, finalPoint3)
	print 'between 2 and 3: ', euclideanDist(finalPoint3, finalPoint2)

	plot(x, y, z)


def plot(x, y, z):
	# Plots path of points
	# Restructure so array values assigned to every N point
	ptx1 = x[0::3]
	pty1 = y[0::3]
	ptz1 = z[0::3]
	ptx2 = x[1::3]
	pty2 = y[1::3]
	ptz2 = z[1::3]
	ptx3 = x[2::3]
	pty3 = y[2::3]
	ptz3 = z[2::3]

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(ptx1, pty1, ptz1, c='red')
	ax.scatter(ptx2, pty2, ptz2, c='blue')
	ax.scatter(ptx3, pty3, ptz3, c='green')

	plt.show()
	

def putInPymol(array):
	# Organizes array into Pymol readable format

	# Get dimensions of array and flatten to make life easier
	x = np.array(array)
	flatArray = x.flatten()
	size = flatArray.shape

	# Round floats
	for n in range(size[0]):
		flatArray[n] = round(flatArray[n], 3)

	# Cast array as strings
	stringArray = map(str, flatArray)

	# Want format: ATOM      0  H0   U 0  10      x y z
	# Must iterate every 3 as we've flattened the array
	file = open('testfile.pdb', 'w')
	for n in range(0, size[0] - 2, 3):
		file.write ('ATOM      ' + str(n / 3) + '  H' + str(n / 3) 
					+ '   U 0  10      ' )
		file.write(stringArray[n] + " "
				 + stringArray[n + 1] + " "
				 + stringArray[n + 2])
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

	"""
	# Test to see if generateAtoms generates
	# the right size array
	def test_generateAtoms(self):
		for [test_num, expected_size] in self.data1:
		    (actual_array, actual_size) = generateAtoms(test_num)
		    self.assertEqual(actual_size, expected_size)
	"""
	"""
	# Test to see if calculateV generates values
	# Compare to calculated values
	def test_calculateV(self):
		for array in self.easy:
			calculated_V = calculateV(array)
			actual_V = [-7.1181416371322853e-19, 
			-7.1181416371322853e-19, -1.3534136350801918e-19, 
			-1.3534136350801918e-19]
			self.assertEqual(calculated_V, actual_V)
	"""
	"""
	# Test putting arrays into Pymol coords
	def test_putInPymol(self):
   		for array in self.easy:
   			putInPymol(array)
   	"""


   	def test_doItAll(self):
   		doItAll(num)


if __name__ == "__main__":
	unittest.main()