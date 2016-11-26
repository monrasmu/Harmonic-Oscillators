# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum bond distance for n H2 atoms
# Currently only working for 3 atoms

import unittest
import time
import sys
import operator
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

# Variable to be changed as needed
# Number of atoms
# Number of steps (5000 appears to be optimal)
num = 4
steps = 5000


def doItAll(num):
	# It does it all!
	# Generates random atoms and moves them to find mimumum potential
	# Updates a progress bar
	# Then plots points and puts final in PDB format
	# But wait, there's more!
	# Times the program too

	tic = time.clock()

	progress = 0
        for step in range(steps):
            progress = updateProgressBar(progress, step, steps)

	print '\n'

	array = generateAtoms(num)
	arrayPoints, vArray = moveMolecule(array, numSteps=steps)
	plotMolecule(arrayPoints)

	print min(vArray), " kJ/mol is the minimum potential achieved"

	toc = time.clock()
	print (toc - tic), "seconds to run full function"

def updateProgressBar(progress, step, steps):
    # Progress bar
    # Taken from EECS183 final project 
    progress += 1
    hashes = ("#" * ((progress * 50) / steps)).ljust(50)
    progressBar = "\r[ " + hashes + " ] Operation in progress..."
    sys.stdout.write(progressBar)
    sys.stdout.flush()
    return progress

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
		index = tree.query_ball_point(x=n, r=2)
		for i in index:
			radius.append(euclideanDist(tree.data[i], n))
	
	return radius


def functionV(r):
		return (4 * _e) * (((_d / r) ** 12) - ((_d / r) ** 6))


def sumV(array):
	# Uses NN search to calculate sum of potential of system
	V = []

	radius  = NNSearch(array)
	
	# Calculate potential for NNs
	for r in radius:
		V.append(functionV(r))
	return sum(V)


"""
def eulerMethod(function, array, numsteps):
	stepSize = 0.1
	v0 = function(array)

	vArray = []
	arrayPoints = []

	# Pass an array of points
	# Function returns potential of system and the array of points
	for n in range(1, numsteps):
		m = function(array)
		v1 = v0  + stepSize * m
		array1 = [(x + stepSize) for x in array]
		arrayPoints.append(array)
		vArray.append(v1)
		array = array1
		v0 = v1
	return vArray, arrayPoints
"""


def moveMolecule(array, numSteps):
	# Moves molecules in random directions by adding random array
	# Will not move points if sum of potential increases
	points = []
	vArray = []

	for i in range(numSteps):
		addArray = np.random.uniform(low=-0.1, high=0.1, size=(len(array),3))
		if sumV(np.add(array, addArray)) < sumV(array):
			array = np.add(array, addArray)
		else: 
			array = array
		points.append(array)
		vArray.append(sumV(array))

	return points, vArray


def plotMolecule(array):
	# Setup for plotting points
	# Puts final output into Pymol and calls plot to plot

	points = []
	finalPointN = []
	x = []
	y = []
	z = []

	# Extract points for plotting
	points = np.squeeze(array)
	for i in points:
		for j in i:
			x.append(j[0])
			y.append(j[1])
			z.append(j[2])

	for n in range(num):
		finalPointN.append((x[-n], y[-n], z[-n]))	

	putInPymol(finalPointN)

	for n in range(num):
		for m in range(1, num):
			if n != m:
				print 'between %i and %i: '  % ((n + 1), (m + 1)), euclideanDist(finalPointN[n], finalPointN[m])

	#plot(x, y, z)
	

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
	# NEED TO FIX TO ACCEPT NEGATIVES

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

"""
class Test(unittest.TestCase):
	data = [([0, 0, 0],
			 [1, 1, 0],
			 [3, 1, 0])]


   	def test_doItAll(self):
   		doItAll(num)
"""

if __name__ == "__main__":
	doItAll(num)
	#unittest.main()