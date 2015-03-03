def makeMultiFrameSimFile(fname, numFrames):
	'''Read in multiple files with the name fname + number + .dat.  Save all those frames to a single file with one header. '''

	import numpy as np
	import struct

	newFile = open(fname + 'allFrames.dat', 'wb')
	#Read in the first file and use this one to write the header.

	#Read in old file
	fIn = open(fname + str(0) + '.dat', 'rb')
	freqstep =float( np.fromfile(fIn, np.double,1) )
	points = int( np.fromfile(fIn, np.int32,1) )
	lines = int( np.fromfile(fIn, np.int32,1) )
	tempReal = np.fromfile(fIn,np.double, points*lines )
	tempImag = np.fromfile(fIn, np.double, points*lines )
	fIn.close()

	#now write the header in the new file
	header = struct.pack('diii', freqstep, points, lines, int(numFrames))
	newFile.write(header)	

	#also write the real then the imaginary part of the file into
	#the header with the frequency index increasing the fastest
	for l in range(lines):
		for p in range(points):
			newFile.write(struct.pack('d', float(tempReal[p + l*points])))

	
	for l in range(lines):
		for p in range(points):
			newFile.write(struct.pack('d', float(tempImag[p + l*points])))
	
	
	for f in range(1,numFrames):

		#Read in old file
		fIn = open(fname + str(f) + '.dat', 'rb')
		freqstep =float( np.fromfile(fIn, np.double,1) )
		points = int( np.fromfile(fIn, np.int32,1) )
		lines = int( np.fromfile(fIn, np.int32,1) )
		tempReal = np.fromfile(fIn,np.double, points*lines )
		tempImag = np.fromfile(fIn, np.double, points*lines )
		fIn.close()

		#write the real then the imaginary part of the file into
		#the file with the frequency index increasing the fastest
		for l in range(lines):
			for p in range(points):
				newFile.write(struct.pack('d', float(tempReal[p + l*points])))

		
		for l in range(lines):
			for p in range(points):
				newFile.write(struct.pack('d', float(tempImag[p + l*points])))
