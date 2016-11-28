"# Harmonic-Oscillators" 

Project studying the Lennard-Jones potential of N atoms.

Requires: unittest, numpy, scipy, matplotlib, mpl_toolkits.mplot3d

Current status: Working for N atoms. Updates progress bar to show time to run. Animation is pretty slow; found it would be simpler to plot initial and final points. Also, now plots potential between two points. Still need to fix PutInPymol

Thoughts: There are definitely better ways to run such minimizations--gradient, etc. This method appears to give a successful "quick and dirty" approach--was fairly quick to code. I need to improve upon memory and time. I initialize a lot of arrays that do not need to be initialized. As well, I should work on splitting up programs into calculating and plotting. Organization could be greatly improved.