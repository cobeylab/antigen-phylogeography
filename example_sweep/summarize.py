#!/usr/bin/env python

import os
import sys
import csv
import sqlite3
import json
import numpy as np
import subprocess
from collections import OrderedDict

def summarize(resultsDb, parName1, parName2, parName3, parName4, parName5, parValue1, parValue2, parValue3, parValue4, parValue5):
	runIds = [row[0] for row in resultsDb.execute(
		"SELECT runId FROM parameters WHERE {0} = ? AND {1} = ? AND {2} = ? AND {3} = ? AND {4} = ?".format(parName1, parName2, parName3, parName4, parName5),
		[parValue1, parValue2, parValue3, parValue4, parValue5]
	)]
	extinctions = list()
	excessDiversities = list()
	fluLike = list()
	
	sideMuRate = list()
	sideMuSize = list()
	sideAgFlux = list()
	trunkMuRate = list()
	trunkMuSize = list()
	trunkAgFlux = list()	
	meanFluxRate = list()
	
	northTrunkPro = list()
	tropicsTrunkPro = list()
	southTrunkPro = list()
	
	northAgLag = list()
	tropicsAgLag = list()
	southAgLag = list()

	
	for runId in runIds:
		try:
			extinctions.append(float(resultsDb.execute('SELECT extinct FROM results WHERE runId = ?', [runId]).next()[0]))
			excessDiversities.append(float(resultsDb.execute('SELECT excessDiversity FROM results WHERE runId = ?', [runId]).next()[0]))
			fluLike.append(float(resultsDb.execute('SELECT fluLike FROM results WHERE runId = ?', [runId]).next()[0]))
			
			sideMuRate.append(float(resultsDb.execute('SELECT sideMuRate FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			sideMuSize.append(float(resultsDb.execute('SELECT sideMuSize FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			sideAgFlux.append(float(resultsDb.execute('SELECT sideAgFlux FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			trunkMuRate.append(float(resultsDb.execute('SELECT trunkMuRate FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			trunkMuSize.append(float(resultsDb.execute('SELECT trunkMuSize FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			trunkAgFlux.append(float(resultsDb.execute('SELECT trunkAgFlux FROM mutationSizeFlux WHERE runId = ?', [runId]).next()[0]))
			
			meanFluxRate.append(float(resultsDb.execute('SELECT meanFluxRate FROM meanFluxRate WHERE runId = ?', [runId]).next()[0]))
			
			northTrunkPro.append(float(resultsDb.execute('SELECT north FROM trunkProportions WHERE runId = ?', [runId]).next()[0]))
			tropicsTrunkPro.append(float(resultsDb.execute('SELECT tropics FROM trunkProportions WHERE runId = ?', [runId]).next()[0]))
			southTrunkPro.append(float(resultsDb.execute('SELECT south FROM trunkProportions WHERE runId = ?', [runId]).next()[0]))

			northAgLag.append(float(resultsDb.execute('SELECT north FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0]))
			tropicsAgLag.append(float(resultsDb.execute('SELECT tropics FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0]))
			southAgLag.append(float(resultsDb.execute('SELECT south FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0]))

			
		except:
			pass
	
	resultsDb.execute(
		'INSERT INTO sweep_summary VALUES (?,?,?,?,?,?,?,?)',
		[parValue1, parValue2, parValue3, parValue4, parValue5, np.sum(extinctions), np.sum(excessDiversities), np.sum(fluLike)]
	)
	
	resultsDb.execute(
		'INSERT INTO agSummaryMedians VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
		[parValue1, parValue2, parValue3, parValue4, parValue5, 
		np.median(sideMuRate), np.median(sideMuSize), np.median(sideAgFlux),
		np.median(trunkMuRate), np.median(trunkMuSize), np.median(trunkAgFlux), 
		np.median(meanFluxRate), 
		np.median(northTrunkPro), np.median(tropicsTrunkPro), np.median(southTrunkPro), 
		np.median(northAgLag), np.median(tropicsAgLag), np.median(southAgLag)]
	)
	
	resultsDb.execute(
		'INSERT INTO agSummaryMean VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
		[parValue1, parValue2, parValue3, parValue4, parValue5, 
		np.mean(sideMuRate), np.mean(sideMuSize), np.mean(sideAgFlux),
		np.mean(trunkMuRate), np.mean(trunkMuSize), np.mean(trunkAgFlux), 
		np.mean(meanFluxRate), 
		np.mean(northTrunkPro), np.mean(tropicsTrunkPro), np.mean(southTrunkPro), 
		np.mean(northAgLag), np.mean(tropicsAgLag), np.mean(southAgLag)]
	)
	
	resultsDb.execute(
		'INSERT INTO agSummarySD VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
		[parValue1, parValue2, parValue3, parValue4, parValue5, np.std(sideMuRate), np.std(sideMuSize), np.std(sideAgFlux),
		np.std(trunkMuRate), np.std(trunkMuSize), np.std(trunkAgFlux), 
		np.std(meanFluxRate), 
		np.std(northTrunkPro), np.std(tropicsTrunkPro), np.std(southTrunkPro), 
		np.std(northAgLag), np.std(tropicsAgLag), np.std(southAgLag)]
	)
	
	resultsDb.execute(
		'INSERT INTO AgLagTrunkPro VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
		[parValue1, parValue2, parValue3, parValue4, parValue5, 
		np.mean(tropicsTrunkPro), np.std(tropicsTrunkPro), 
		np.mean(tropicsAgLag), np.std(tropicsAgLag),
		np.sum(extinctions), np.sum(excessDiversities), np.sum(fluLike)]
	)
	
	
	
if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	
	parName1 = "relativeN"
	parName2 = "tropicFractionI0"
	parName3 = "relativeR0"
	parName4 = "relativeTurnover"
	parName5 = "seasonalAmplitude"
	
	resultsDb = sqlite3.connect('results.sqlite')
	
	parameterSets = [row for row in resultsDb.execute(
		"SELECT DISTINCT {0}, {1}, {2}, {3}, {4} FROM parameters".format(parName1,parName2,parName3,parName4,parName5)
	)]
	
	# Create table containing various metrics for each rate, sd pair
	resultsDb.execute(
		"CREATE TABLE sweep_summary ({0},{1},{2},{3},{4},{5},{6},{7})".format(parName1,parName2,parName3,parName4,parName5,"extinctions","excessDiversities","fluLike")
	)
	
	# Create table for antigenic summary
	resultsDb.execute(
		"CREATE TABLE agSummaryMedians ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17})".format(parName1,parName2,parName3,parName4,parName5,"sideMuRate","sideMuSize","sideAgFlux","trunkMuRate","trunkMuSize","trunkAgFlux","meanFluxRate","northTrunkPro","tropicsTrunkPro","southTrunkPro","northAgLag","tropicsAgLag","southAgLag")
	)
	
	resultsDb.execute(
		"CREATE TABLE agSummaryMean ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17})".format(parName1,parName2,parName3,parName4,parName5,"sideMuRate","sideMuSize","sideAgFlux","trunkMuRate","trunkMuSize","trunkAgFlux","meanFluxRate","northTrunkPro","tropicsTrunkPro","southTrunkPro","northAgLag","tropicsAgLag","southAgLag")
	)
	
	resultsDb.execute(
		"CREATE TABLE agSummarySD ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17})".format(parName1,parName2,parName3,parName4,parName5,"sideMuRate","sideMuSize","sideAgFlux","trunkMuRate","trunkMuSize","trunkAgFlux","meanFluxRate","northTrunkPro","tropicsTrunkPro","southTrunkPro","northAgLag","tropicsAgLag","southAgLag")
	)
	
	resultsDb.execute(
		"CREATE TABLE AgLagTrunkPro ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11})".format(parName1,parName2,parName3,parName4,parName5,"tropicsTrunkPro","tropicsTrunkProSD","tropicsAgLag","tropicsAgLagSD","extinctions","excessDiversities","fluLike")
	)

	for parValue1, parValue2, parValue3, parValue4, parValue5 in parameterSets:
		summarize(resultsDb, parName1, parName2, parName3, parName4, parName5, parValue1, parValue2, parValue3, parValue4, parValue5)
	
	resultsDb.commit()
	
