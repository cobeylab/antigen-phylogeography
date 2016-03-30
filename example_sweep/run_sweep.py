#!/usr/bin/env python

# jobs are indexed staring from 0

import csv
import os
import json
import numpy as np
import subprocess
from collections import OrderedDict

###############################################################################

### Modify these values and the generateJobs() function ###

# Set DRY = False to actually submit jobs
DRY = False

maxJobs = 500
steps = 12

R0S = np.linspace(1.2,2.8,steps)

NUM_RUNS = 20
totalJobs = NUM_RUNS*len(R0S)
if totalJobs < maxJobs:
	maxJobs = totalJobs
N_PER_JOB = totalJobs/maxJobs+1

totalN = 15000000
birthRate = 0.000091
Reff = 1.0
R0 = 1.8
nu = 0.2
beta = R0 * (nu+birthRate)
initialTraitA = -4.0
smithConversion = 0.07
muPhenotype = 0.0001
betweenDemePro = 0.001
deltaT = 0.1
startAtEquilibriumInfected = [False]
startAtEquilibriumImmune = [False]
meanStep = 0.6
sdStep = 0.3


demeCount = 1
demeNames = ["pop"]
initialNs = [15000000]
demeOffsets = [0]

relativeN = 1.0
tropicFractionI0 = 1.0
relativeR0 = 1.0
relativeTurnover = 1.0
seasonalAmplitude = 0.0

def generateJobs():	
	jobNum = 0

	for R0 in R0S:
		
		
		birthRates = [birthRate]
		deathRates = birthRates
		
		beta = R0 * (nu + birthRates[0])

		demeAmplitudes = [0.0]
		
		demeBaselines = [1.0]
		
		#start at equilibrium
	
		initialPrS = Reff/(R0)

		initialPrI = birthRates[0]/beta*((R0)-1.0)
		totalInitialIs = int(initialPrI * totalN)

		initialPrR = (1.0 - initialPrS - initialPrI)/(1.0-smithConversion*abs(initialTraitA))
		initialIs = [totalInitialIs]
		initialNs = [totalN]
	
		for runNum in range(NUM_RUNS):
			# This is the name that SLURM uses to identify the job
			jobName = jobNum

			# This is the subdirectory inside 'results' used to run the job
			jobSubdir = '{}'.format(jobName)

			paramDict = OrderedDict([
				('smithConversion',smithConversion),
				('meanStep',meanStep),
				('sdStep',sdStep),
				('initialIs',initialIs),
				('initialPrR',initialPrR),
				('totalN',totalN),
				('birthRate',birthRates),
				('deathRate',deathRates),
				('muPhenotype',muPhenotype),
				('seasonalAmplitude', seasonalAmplitude),
				('demeAmplitudes', demeAmplitudes),
				('demeBaselines', demeBaselines),
				('initialNs',initialNs),
				('relativeN', relativeN),
				('tropicFractionI0', tropicFractionI0),
				('relativeR0', relativeR0),
				('relativeTurnover', relativeTurnover),
				('betweenDemePro',betweenDemePro),
				('deltaT',deltaT),
				('startAtEquilibriumImmune',startAtEquilibriumImmune),
				('startAtEquilibriumInfected',startAtEquilibriumInfected),
				('R0',R0),
				('initialTraitA',initialTraitA),
				('beta',beta),
				('demeNames',demeNames),
				('demeOffsets',demeOffsets),
				('demeCount',demeCount),
			])
			yield (jobName, jobSubdir, paramDict)
			jobNum += 1


###############################################################################

def writeParameters(resultsDir, jobSubdir, paramDict):
	jobDir = os.path.join(resultsDir, jobSubdir)
	os.makedirs(jobDir)
	paramsFilename = os.path.join(jobDir, 'parameters.json')
	with open(paramsFilename, 'w') as paramsFile:
		json.dump(paramDict, paramsFile, indent=2)
		paramsFile.write('\n')

def submitJob(rootDir, sweepDir, resultsDir):
	# Construct SLURM command
	sbatchFilename = os.path.join(sweepDir, 'job.sbatch')
	submitCommand = [
		'sbatch',
		'-J{0}'.format('Antigen'),
		'-D{0}'.format(resultsDir),
		'--array=0-{0}'.format(maxJobs-1),
		sbatchFilename
	]
	
	# Print command to terminal
	process = subprocess.Popen(['echo'] + submitCommand)
	process.wait()
	
	# Construct environment variables
	env = dict(os.environ)
	env['ANTIGEN_ROOT'] = rootDir
	env['N_PER_JOB'] = str(N_PER_JOB)
	
	# Actually run command
	if not DRY:
		process = subprocess.Popen(submitCommand, env=env)
		process.wait()

def loadUncommentedJsonString(filename):
	lines = list()
	with open(filename) as jsonFile:
		for line in jsonFile:
			lines.append(line.rstrip('\n').split('//')[0])
	return '\n'.join(lines)

if __name__ == '__main__':
	sweepDir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
	rootDir = os.path.abspath(os.path.join(sweepDir, '..'))
	resultsDir = os.path.join(sweepDir, 'results')
	
	cParamsFilename = os.path.join(sweepDir, 'constant_parameters.json')
	cParamsJson = loadUncommentedJsonString(cParamsFilename)
	cParamsDict = json.loads(cParamsJson, object_pairs_hook=OrderedDict)
	
	for jobName, jobSubdir, paramDict in generateJobs():
		allParamsDict = OrderedDict()
		allParamsDict.update(cParamsDict)
		allParamsDict.update(paramDict)
		writeParameters(resultsDir, jobSubdir, allParamsDict)
		
	submitJob(rootDir, sweepDir, resultsDir)
