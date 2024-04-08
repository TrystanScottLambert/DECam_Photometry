##################################
#
# Script to make DECAM script 
#
##################################

import numpy as np 
from astropy.coordinates import SkyCoord 
import astropy.units as u
import json
import datetime

def prettifySkyCoords(SkyCoordObject):
	stringWithLetters = SkyCoordObject.to_string('hmsdms')
	letters = ['h','m','d']
	for letter in letters:
		stringWithLetters = stringWithLetters.replace(letter,':')

	stringWithNoLetters = stringWithLetters.replace('s','')
	return stringWithNoLetters

class ExposureSequence():
	def __init__(self,ra,dec,filterChoice,exptime,numberOfExposures):
		self.ra = (int(ra.split(':')[0]) + int(ra.split(':')[1])/60 + float(ra.split(':')[2])/3600)*(360/24)
		
		if dec[0] == '-':
			self.dec = int(dec.split(':')[0]) - int(dec.split(':')[1])/60 - float(dec.split(':')[2])/3600
		else:
			self.dec = int(dec.split(':')[0]) + int(dec.split(':')[1])/60 + float(dec.split(':')[2])/3600
		
		self.expTime = exptime
		self.filterChoice = filterChoice
		self.numberExposures = numberOfExposures

	def dither(self):
		decOffset = 5/60  # 5 arminute offset to play in center of field
		ditherOffset = 120/3600 # 120 arcsecond dithers +-
		raDithers = np.random.random(self.numberExposures)*2*ditherOffset - ditherOffset
		decDithers = np.random.random(self.numberExposures)*2*ditherOffset - ditherOffset
		self.raPositions = self.ra + raDithers 
		self.decPositions = self.dec + decDithers + decOffset

	def addSequence(self,fileVariable,includeEndComma=True):
		self.dither()
		for i in range(len(self.raPositions)):
			c = SkyCoord(ra=self.raPositions[i]*u.deg, dec=self.decPositions[i]*u.deg)
			cString = prettifySkyCoords(c)
			raToWrite = cString.split(' ')[0]
			decToWrite = cString.split(' ')[1]
			fileVariable.write(' { \n')
			fileVariable.write('  "expType": "object", \n')
			fileVariable.write('  "object": "VIK_J2348-3054", \n')
			fileVariable.write(f'  "RA": "{raToWrite}", \n')
			fileVariable.write(f'  "dec": "{decToWrite}", \n')
			fileVariable.write(f'  "filter": "{self.filterChoice}", \n')
			fileVariable.write(f'  "expTime": {self.expTime} \n')
			if i == self.numberExposures-1 and includeEndComma == False:
				fileVariable.write(' } \n')
			else:
				fileVariable.write(' }, \n')

def calculateTotalTime(scriptName):
	f = open(scriptName)
	data = json.load(f)
	times = [val['expTime'] for val in data]
	return np.sum(times) + 28*len(times)

def printTotalTime(scriptName):
	seconds = int(calculateTotalTime(scriptName))
	print(str(datetime.timedelta(seconds=seconds)))
	
def constructScript(filename,*args):
	f = open(filename,'w')
	f.write('[ \n')
	for i in range(len(args)):
		if i == len(args)-1:
			args[i].addSequence(f,includeEndComma=False)
		else:
			args[i].addSequence(f)
	f.write('] \n')
	f.close()

	print('Total Time of Script ~ ', str(datetime.timedelta(seconds=int(calculateTotalTime(filename)))))

if __name__ == '__main__':
	RA =  '23:48:33.34'
	Dec = '-30:54:10.0'
	NarrowSequence_1  = ExposureSequence(RA,Dec,'i',200,10)
	#NarrowSequence_1  = ExposureSequence(RA,Dec,'i',200,10)
	constructScript('AssefLambert_i_band7.json',NarrowSequence_1)#,BroadSequence_1,NarrowSequence_2,BroadSequence_2)


23:48:33.34 -30:54:10.0