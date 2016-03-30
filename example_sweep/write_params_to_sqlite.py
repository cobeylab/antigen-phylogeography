#!/usr/bin/env python

import os
import sys
import csv
import sqlite3
import json
import numpy as np
import subprocess
from collections import OrderedDict

def insert(db, table, values):
    return db.execute(
        '''INSERT INTO {0} VALUES ({1})'''.format(
            table,
            ','.join(['?'] * len(values))
        ),
        values
    )
    
def loadJob(comboDb, resultsDir, runId, jobDir):
	# Load dictionary of parameter values from parameters.json file
	jsonFile = os.path.join(jobDir, 'parameters_out.json')
	if(os.path.isfile(jsonFile)):
		paramsJson = loadUncommentedJsonString(jsonFile)
		paramsDict = json.loads(paramsJson, object_pairs_hook=OrderedDict)
		paramVals = paramsDict
	
		insert(comboDb, 'parameters', [runId] + [str(paramVals[paramName]) for paramName in paramsDict.keys()])
		comboDb.commit()
	
def loadJobs(comboDb, sweepDir, resultsDir):
	# Extract parameters names
	subDirs = os.listdir(resultsDir)
	runId = 0
	for subDir in subDirs:
		print subDir
		jobDir = os.path.join(resultsDir, subDir)
		runId = subDir
		if os.path.isdir(jobDir):
			loadJob(comboDb, resultsDir, runId, jobDir)
			
def loadUncommentedJsonString(filename):
	lines = list()
	with open(filename) as jsonFile:
		for line in jsonFile:
			lines.append(line.rstrip('\n').split('//')[0])
	return '\n'.join(lines)

if __name__ == "__main__":
	sweepDir = os.path.abspath(os.path.dirname(__file__))
	resultsDir = os.path.join(sweepDir, 'results')
	
	paramsJson = loadUncommentedJsonString(os.path.join(resultsDir,'0/parameters_out.json'))
	paramsDict = json.loads(paramsJson, object_pairs_hook=OrderedDict)
	paramVals = paramsDict
	
	comboDb = sqlite3.connect(os.path.join(sweepDir, 'parameters_out.sqlite'))
	
	# Index of runs with parameters
	comboDb.execute(
		"CREATE TABLE parameters (runId INTEGER, {0})".format(
			', '.join([paramName + ' ' + 'TEXT' for paramName in paramsDict.keys()])
		)
	)
	
	comboDb.commit()
	
	loadJobs(comboDb, sweepDir, resultsDir)

	comboDb.close()
	