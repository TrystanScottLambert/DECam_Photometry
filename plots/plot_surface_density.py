""" making the onsky surface density plot """

import numpy as np 
import pylab as plt 

INFILE = 'candidates.txt'

ra_candidates, dec_candidates = np.loadtxt(INFILE,unpack=True)