'''
Make Oligo sequences for human TF library
Input: database of TFs, after finding transcripts and then extracting a CDS for each gene
Output: oligo sequences

Josh Tycko
October 2018

'''

import pandas as pd
import numpy as np
import math
from itertools import islice
import argparse
import sys
from dnachisel import *

parser = argparse.ArgumentParser(description='Get sequence of best transcript')
parser.add_argument('in_file', type=str,
                    help='Input file')
parser.add_argument('out_file', type=str,
                    help='Output filename start')

# Optional arguments:
parser.add_argument('-o', '--optimize', dest='opt',
                    help='Optimize the DNA sequences',
                    type=str, default='')
parser.add_argument('-s', dest = 'step', type=str,
                    help='sliding window step size', default = 'scan')

'''
Scan mode: can try many step sizes, do not optimize
Optimize mode: slower, so best to do this only for one selected step size
'''

def optimizeOligo(dna_sequence, pattern):
	problem = DnaOptimizationProblem(
			sequence=dna_sequence,
			constraints=[AvoidPattern(pattern),
			             EnforceGCContent(mini=0.20, maxi=0.75, window=50),
			             EnforceTranslation()],
			objectives=[CodonOptimize(species='h_sapiens')]
			)
	try:
		problem.resolve_constraints()
	except (TypeError, KeyError):
		print('optimization failed')
		return dna_sequence
	problem.optimize()
	optDNA = str(problem.sequence)
	return optDNA

def slidingWindow(sequence,winSize,step=1):
	"""Returns a generator that will iterate through
	the defined chunks of input sequence.  Input sequence
	must be iterable."""

	sequence = sequence[:len(sequence)]

	# Verify the inputs
	try: it = iter(sequence)
	except TypeError:
	    raise Exception("**ERROR** sequence must be iterable.")
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
	    raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
	    raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > len(sequence):
	    raise Exception("**ERROR** winSize must not be larger than sequence length.")

	# Pre-compute number of chunks to emit
	numOfChunks = ((len(sequence)-winSize)/step)+1

	# Do the work
	for i in range(0,numOfChunks*step,step):
		yield sequence[i:i+winSize]
	# Add one more chunk, which goes to the very end
	yield sequence[-winSize:] 

args = parser.parse_args()
TFdf = pd.read_csv(args.in_file)

oligoLen = 300 - 2*(11+19) # 11 for TypeIIs restriction cloning and 19 for primer
print(oligoLen/3, 'amino acids per oligo')
# step = args.step
# if step%3 != 0:
# 	print('step size must be in-frame mod3')
# 	sys.exit()

def makeOligos(df, step, filterDesc):
	rows_list = []
	for index, row in df.iterrows():
		if pd.isna(row['cds']): # Skip the rows that have no CDS
			continue
		row_dict = {}
		cds = row['cds']
		cds = cds[:-3]  # don't include STOP codon
		# numChunks = math.ceil(len(cds)/oligoLen)
		if len(cds) <= oligoLen:
			row_dict.update({'Ensembl ID': row['Ensembl ID'],'Segment': 0,'Sequence': cds})
			rows_list.append(row_dict)
			continue

		chunks = ["".join(x) for x in slidingWindow(cds, winSize = oligoLen, step = step)]
		segment = 0
		for chunk in chunks:
			row_dict.update({'Ensembl ID': row['Ensembl ID'],'Segment': segment,'Sequence': chunk})
			rows_list.append(row_dict)
			row_dict = {}
			segment += 1
	oligodf = pd.DataFrame(rows_list, columns = ['Ensembl ID', 'Segment', 'Sequence'])

	## Remove duplicates:
	dups = sum(oligodf['Sequence'].duplicated())
	if dups > 0:
		oligodf['Ensembl ID with redundant sequence'] = None
		for cds in set(list(oligodf[oligodf['Sequence'].duplicated(keep=False)]['Sequence'])): # this makes it way faster
			EnsemblIDs = list(oligodf.loc[oligodf['Sequence'] == cds]['Ensembl ID'])
			if EnsemblIDs[0] == EnsemblIDs[1]:
				continue
			if len(EnsemblIDs) > 1:
				oligodf.loc[oligodf['Sequence'] == cds, 'Ensembl ID with redundant sequence'] = str(EnsemblIDs)

		oligodf = oligodf.drop_duplicates('Sequence')

	## Optimize to remove BsmBI and fix GC content
	if args.opt == 'CodonOpt':
		BsmBIpattern = enzyme_pattern("BsmBI")
		for cds in set(list(oligodf['Sequence'])):
			try:
				optcds = optimizeOligo(cds, BsmBIpattern)
			except KeyError:
				continue
			oligodf.loc[oligodf['Sequence'] == cds, 'Sequence'] = optcds


	print(filterDesc, str(step), oligodf['Segment'].count(), 'Duplicates', dups)
	oligodf.to_csv(args.out_file + filterDesc +'_step' + str(step) + args.opt + '.csv', index = False)

steplist = []
if args.step == 'scan':
	steplist = [180,150,120,90,60,30,15]
else:
	chosenStep = int(args.step)
	steplist.append(chosenStep)

# filterDesc = 'all'
# for step in steplist:
# 	makeOligos(TFdf, step, filterDesc)

# curdf = TFdf.loc[TFdf['Is TF?'] == 'Yes']
# filterDesc = 'Lambert et al. curated'
# for step in steplist:
# 	makeOligos(curdf, step, filterDesc)

joshdf = TFdf.loc[TFdf['Josh subset'] == 'Yes']
filterDesc = 'Josh subset'
for step in steplist:
	makeOligos(joshdf, step, filterDesc)
