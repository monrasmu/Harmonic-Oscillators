# Harmonic Oscillator N-D
# Uses LJ potential to calculate minimum
# bond distance for n H2 atoms

import unittest
import cmath
import numpy as np

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
	radius = cmath.sqrt((array[1][1] - array[2][1]) ** 2 
					  + (array[1][2] - array[2][2]) ** 2 
					  + (array[1][3] - array[2][3]) ** 2)
	_V = ((4 * _e) * (((_d / radius) ** 12) - ((_d / radius) ** 6)))
	return _V

def putInPymol(array):
	print array
	# First cast array as strings
	str_arr = [str(num) for num in array]
	
	# put into PyMol readable format
	file = open('testfile.pdb', 'w')
	file.write ('ATOM      1  H1   U 0  10      ' 
				+ strx1 + '.000  ' + stry1 + '.000  ' 
				+ strz1 + '.000 \n'
			  + 'ATOM      2  H2   U 0  10      ' 
				+ strx2 + '.000  ' + stry2 + '.000  ' 
				+ strz2 + '.000 \n'
			  + 'ATOM      3  H3   U 0  10      ' 
				+ strx3 + '.000  ' + stry3 + '.000  ' 
				+ strz3 + '.000 \n')
	file.close()



class Test(unittest.TestCase):
    # Test cases
    data1 = [(2, (2, 3)),
    		(3, (3, 3)),
    		(0, (0, 3))]
    data2 = 
    data3 = [[2, 3],
    		 [3, 4, 4]]

    # test to see if generateAtoms generates
    # the right size array
    def test_generateAtoms(self):
        for [test_num, expected_size] in self.data1:
            (actual_array, actual_size) = generateAtoms(test_num)
            self.assertEqual(actual_size, expected_size)

    # test to see if calculateV generates values
    # COMPARE TO CALCULATED?
    def test_calculateV(self):
   		for array in self.data2:

   	def test_putInPymol(self):
   		for array in self.data3:
   			putInPymol(array)




if __name__ == "__main__":
	unittest.main()