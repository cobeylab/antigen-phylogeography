#!/usr/bin/env python

import os
import sys
import csv
import json
import subprocess
import shutil

def loadUncommentedJsonString(filename):
	lines = list()
	with open(filename) as jsonFile:
		for line in jsonFile:
			lines.append(line.rstrip('\n').split('//')[0])
	return '\n'.join(lines)

if __name__ == '__main__':
	withRestarts = False
	if len(sys.argv) > 1:
		if sys.argv[1] == '--with-restarts':
			withRestarts = True
	
	rootDir = os.path.dirname(__file__)
	
	args = [
		'java',
		'-ea',
		'-Xmx24G',
		'-XX:+UseSerialGC',
		'-cp',
		os.path.join(rootDir, 'lib', '*'),
		'antigen.Antigen'
	]
	
	params = json.loads(loadUncommentedJsonString('parameters.json'))
	endDay = float(params['endDay'])
	printStep = float(params['printStep'])
	
	done = False
	while not done:
		process = subprocess.Popen(args)
		process.wait()
		
		if withRestarts:
			restartLim = 10.0
			restartCount = 0.0
			with open('out.timeseries') as f:
				cr = csv.reader(f, delimiter='\t')
				for row in cr:
					pass
				lastDay = float(row[0]) * 365.0
			
				if (abs(lastDay - endDay) < printStep + 1.0 or restartCount >= restartLim):
					done = True
				else:
					print('Simulation aborted; restarting')
					restartCount = restartCount + 1.0
		else:
			done = True
	print('Simulation completed')
