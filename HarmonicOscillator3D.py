# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum bond distance for n H2 atoms

import unittest
import time
import sys
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
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
# Number of steps
num = 20
steps = 1000


def doItAll(num):
	# It does it all!
	# Generates random atoms and moves them to find mimumum potential
	# Updates a progress bar
	# Then plots points and puts final in PDB format
	# But wait, there's more!
	# Times the program too
	newlist = []
	tic = time.clock()

	array = generateAtoms(num)
	arrayPoints, vArraySum = moveMolecule(array, numSteps=steps)
	plotSetUp(arrayPoints)
	print min(vArraySum)," kJ/mol is the minimum potential achieved for the system"
	toc = time.clock()
	print round((toc - tic), 3), "seconds to run"


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
	array = np.random.uniform(low=-0.001, high=0.001, size=(num,3))
	array_size = array.shape
	array = np.squeeze(array)
	return array


def euclideanDist(point1, point2):
	# Returns euclidean distance between two points
	return spatial.distance.euclidean(point1, point2)


def NNSearchSetUp(array):
	# Takes array and assigns them to tuples
	return tuple(map(tuple, array))


def NNSearch(points):
	# Query tree for nearest neighbors using KD tree
	# Requires that points are in tuples
	# Returns radius as list of NN distances
	radius = []
	for n in points:
		output = [i for i in points if i != n]
		tree = spatial.KDTree(output)
		index = tree.query_ball_point(x=n, r=2)
		for i in index:
			radius.append(euclideanDist(tuple(tree.data[i]), n))
	return radius


def functionV(r):
	# Function to calculate LJ potential
	return ((4 * _e) * (((_d / r) ** 12) - ((_d / r) ** 6)))


def sumV(array):
	# Uses NN search to calculate sum of potential of system
	V = []

	# Performs NN search to collect nearby radii
	points = NNSearchSetUp(array)

	radius  = NNSearch(points)
	
	# Calculate potential for NNs
	for r in radius:
		V.append(functionV(r))
	return sum(V)


def vArrays(array):
	# Uses NN search to calculate the potential of system
	V = []

	# Performs NN search to collect nearby radii
	points = NNSearchSetUp(array)

	radius  = NNSearch(points)
	
	# Calculate potential for NNs
	for r in radius:
		V.append(functionV(r))
	return V

"""
def gradient(function):
	radius = NNSearch(array)
	zs = np.array([functionV(rad) for rad in radius])
	gradx = gradient(zs)
	print gradx
	return np.gradient(function)


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
	# Will not move points if potential of system increases
	points = []
	vArray = []
	vArraySum = []
	progress = 0

	points.append(array)
	for i in range(numSteps):
		addArray = np.random.uniform(low=-0.1, high=0.1, size=(num,3))
		if sumV(np.add(array, addArray)) < sumV(array):
			array = np.add(array, addArray)
		else: 
			array = array

		points.append(array)
		vArray.append(vArrays(array))
		vArraySum.append(sumV(array))

		# Updates progress bar
		progress = updateProgressBar(progress, i, numSteps)

	# Pass atoms to plotV to plot
	plotV(points, vArray)

	return points, vArraySum


def plotV(points, vArray):
	# Plots potential as a function of radius between two points
	radius = []
	vPoints = []
	vArray = np.squeeze(vArray)

	# Gets distance between first two points to plot example V vs radius
	for n in range(steps):
		radius.append(euclideanDist(points[n][0], points[n][1]))
		vPoints.append(vArray[n][0])
	
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	ax.set_title('Potential vs Radius')
	ax.set_xlabel('radius')
	ax.set_ylabel('potential (kJ/mol)')
	ax.scatter(radius, vPoints, s=50)
	ax.set_xlim(0, 2)
	ax.set_ylim(-0.5, 1)

	x = np.linspace(0.001, 2, 1000)
	y = functionV(x)
	ax.plot(x, y)
	

def plotSetUp(array):
	# Setup for plotting points
	# Puts final output into Pymol and calls plot to plot
	# Prints average bond distance !OJO! New fxn
	points = []
	finalPointN = []
	intialPointN = []
	x = []
	y = []
	z = []
	radius = []

	# Extract points for plotting
	points = np.squeeze(array)
	for i in points:
		for j in i:
			x.append(j[0])
			y.append(j[1])
			z.append(j[2])

	for n in range(num):
		finalPointN.append((x[-num + n], y[-num + n], z[-num + n]))
		intialPointN.append((x[n], y[n], z[n]))	
	
	putInPymol(finalPointN)

	"""
	# Prints distance to all neighbors
	for n in range(num):
		for m in range(1, num):
			if n != m:
				# print 'between %i and %i: '  % ((n + 1), (m + 1)), euclideanDist(finalPointN[n], finalPointN[m])
				radii.append(euclideanDist(finalPointN[n], finalPointN[m]))
	"""

	# Prints average bond distance to neighbors
	print '\n'
	print 'Average bond distance is', sum(NNSearch(finalPointN)) / len(NNSearch(finalPointN))

	plot(intialPointN, finalPointN)
	

def plot(initial, final):
	# Plots initial and final output
	# Input must be array of tuples of points
	ptxN = []
	ptyN = []
	ptzN = []
	ptx1 = []
	pty1 = []
	ptz1 = []
	for n in range(num):
		ptxN.append(final[n][0])
		ptyN.append(final[n][1])
		ptzN.append(final[n][2])
		ptx1.append(initial[n][0])
		pty1.append(initial[n][1])
		ptz1.append(initial[n][2])

	fig = plt.figure()
	ax = fig.add_subplot(2, 1, 1, projection='3d')
	ax.set_title('Initial')
	ax.set_xlim3d(-1, 1)
	ax.set_ylim3d(-1, 1)
	ax.set_zlim3d(-1, 1)
	ax2 = fig.add_subplot(2, 1, 2, projection='3d')
	ax2.set_title('Final')
	ax2.set_xlim3d(-1, 1)
	ax2.set_ylim3d(-1, 1)
	ax2.set_zlim3d(-1, 1)
	colors = np.random. rand(num)
	ax.scatter(ptx1, pty1, ptz1, c=colors, s=50)
	ax2.scatter(ptxN, ptyN, ptzN, c=colors, s=50)
	plt.show()


def plotAnimation(x, y, z):
	# Plots path of points
	# Restructure so array values assigned to every 3 points
	# !OJO! Currently only plotting effectively for 3 atoms
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
		ax.scatter(ptx1[:(5*i)], pty1[:(5*i)], ptz1[:(5*i)], s=area, c='red')
		ax.scatter(ptx2[:(5*i)], pty2[:(5*i)], ptz2[:(5*i)], s=area, c='blue')
		ax.scatter(ptx3[:(5*i)], pty3[:(5*i)], ptz3[:(5*i)], s=area, c='green')

	ani = FuncAnimation(fig, animate, frames=len(ptx1), interval=10)
	
	plt.show()
	

def putInPymol(array):
	# Organizes array into Pymol readable format
	# !OJO! Need to fix to accept negative values

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


if __name__ == "__main__":
	doItAll(num)