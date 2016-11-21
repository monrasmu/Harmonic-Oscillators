# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum
# bond distance for n H2 atoms
# ONLY WORKING FOR 3 ATOMS RIGHT NOW

import unittest
import time
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
import matplotlib.lines as line
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

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
	calculateV(array)
	#moveMolecule(array)


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

	def functionV(x, y):
		r = euclideanDist(x, y)
		return (4 * _e) * (((_d / r) ** 12) - ((_d / r) ** 6))

	eulerMethod(functionV, 1, 10)
	return V


def sumV(array):
	return sum(array)

def gradientV(array):
	grad = np.gradient(array)
	print grad
	return grad

def eulerMethod(function, guess, numsteps):
	x1 = guess
	x2 = guess + 0.1
	y1 = function(x1)
	y2 = function(x2)
	dx = (y2 - y1) / (x2 - x1)

	for n in range(numsteps):
		y[n + 1] =y[n] + dx*function(x[n], y[n])
	print y
	return y


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
	timestop = time.time() + 5    # 5 seconds
	while all(V > -0.09 for V in calculateV(array)):
		addArray = np.random.uniform(low=-1, high=1, size=(len(array),3))
		array = np.add(array, addArray)
		points.append(array)
		if time.time() > timestop:
			break

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
	# Restructure so array values assigned to every 3 points
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

	# Increase size of points as steps increase
	area = [8*n for n in range(len(ptx1))]

	# Animate, iterating by 2* i-th element from points
	def animate(i):
		ax.scatter(ptx1[:(2*i)], pty1[:(2*i)], ptz1[:(2*i)], s=area, c='red')
		ax.scatter(ptx2[:(2*i)], pty2[:(2*i)], ptz2[:(2*i)], s=area, c='blue')
		ax.scatter(ptx3[:(2*i)], pty3[:(2*i)], ptz3[:(2*i)], s=area, c='green')

	ani = FuncAnimation(fig, animate, frames=len(ptx1), interval=200)
	
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
	data = [([0, 0, 0],
			 [1, 1, 0],
			 [3, 1, 0])]


   	def test_doItAll(self):
   		doItAll(num)


if __name__ == "__main__":
	unittest.main()