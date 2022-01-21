'''
make deep mutational scanning oligos for any input

Input: a protein sequence
Output: DMS library oligos

Josh Tycko
October 2018, updated Aug 2020
'''

import pandas as pd
# from Bio import SeqIO, Seq
import re
import argparse
from dnachisel import *

parser = argparse.ArgumentParser(description='Make sequences for DMS library')
# parser.add_argument('residues', type=int,
#                     help='Number of residues per element')
parser.add_argument('input_sequence', type=str,
                    help='Input AA sequence')
parser.add_argument('name', type=str,
                    help='Protein name')
parser.add_argument('-s','--swaps', dest='swaps', type=int,
                    help='Max number of consecutive AA swaps (eg. 3 gives single, double, and triple swaps', default=1)
parser.add_argument('-c', '--controls', dest = 'controls', type=int,
                    help='Number of random control sequences to include',
                    default = 500)
parser.add_argument('-p', '--pairs', dest = 'pairwise', type=bool,
                    help='whether to do pairwise mutations for epistasis',
                    default = False)

#########################################################################################################

def optimizeOligo(dna_sequence, pattern):
	dna_sequence = dna_sequence.upper()
	problem = DnaOptimizationProblem(
			sequence=dna_sequence,
			constraints=[AvoidPattern(pattern), AvoidPattern('7xC'),
			             EnforceGCContent(mini=0.20, maxi=0.75, window=40),
			             AvoidRareCodons(species='h_sapiens', min_frequency=0.05),
			             EnforceTranslation()],
			objectives=[CodonOptimize(species='h_sapiens', method='match_codon_usage'), UniquifyAllKmers(15)]
			)
	try:
		problem.resolve_constraints()
	except (TypeError, KeyError):
		print('optimization failed')
		return dna_sequence
	problem.optimize()
	optDNA = str(problem.sequence)
	return optDNA

#########################################################################################################
args = parser.parse_args()
BsmBIpattern = EnzymeSitePattern("BsmBI")
print('Avoiding BsmBI sites')
numControls = args.controls
numWTcontrols = 100
WTseq_AA = args.input_sequence
ProteinName = args.name

rows_list = [] #list of sequences

# First add WT controls
print('Generating up to '+str(numWTcontrols)+' wildtype controls with varied codon usage')
dna_sequence = reverse_translate(WTseq_AA)
i = 0
stop = 0
listWTs = []
while i < numWTcontrols and stop < 1000:
	WTseq_DNA = optimizeOligo(dna_sequence, BsmBIpattern)
	if WTseq_DNA not in listWTs:
		listWTs.append(WTseq_DNA)
		rows_list.append({'label':'WT;'+ProteinName +';' +str(i),'Swap size': 0, 'Variant protein': WTseq_AA, 'Variant DNA': WTseq_DNA, 'Category': 'WT'})
		i += 1
	stop += 1

print('Successfully generated '+str(len(listWTs))+' wildtype controls with varied codon usage')

# Make random negative controls
print('Generating '+str(numControls)+' random controls')
for ctrlNumber in range(0,numControls):
	# get a random protein with no stop codons
	dna_sequence = random_dna_sequence(len(WTseq_AA*3))
	while '*' in translate(dna_sequence):
		dna_sequence = random_dna_sequence(len(WTseq_AA*3))
	# optimize to remove RE sites	
	variantDNA = optimizeOligo(dna_sequence, BsmBIpattern)
	variant = translate(variantDNA)
	rows_list.append({'label':'Random_control;'+str(ctrlNumber), 'Variant protein': variant, 'Variant DNA': variantDNA, 'Category': 'Random control'})

AAdict = {'Alanine':'A',
	'Cysteine':'C',
	'AsparticAcid':'D',
	'GlutamicAcid':'E',
	'Phenylalanine':'F',
	'Glycine':'G',
	'Histidine':'H',
	'Isoleucine':'I',
	'Lysine':'K',
	'Leucine':'L',
	'Methionine':'M',
	'Asparagine':'N',
	'Proline':'P',
	'Glutamine':'Q',
	'Arginine':'R',
	'Serine':'S',
	'Threonine':'T',
	'Valine':'V',
	'Tryptophan':'W',
	'Tyrosine':'Y'}

# Make mutants
print('Generating scan of mutants')
for AAname, AAswap in AAdict.items():
	print('Scanning', AAname)
	for swapNum in range(1, args.swaps+1):
		for pos in range(0, len(WTseq_AA)):
			seq = list(WTseq_AA) # set to original
			if pos+swapNum <= len(seq):
				seq[pos:pos+swapNum] = AAswap*swapNum
				variant = ''.join(seq) # protein sequence
				if variant==WTseq_AA:
					continue
				dna_sequence = reverse_translate(variant)
				variantDNA = optimizeOligo(dna_sequence, BsmBIpattern)
				rows_list.append({'label':'DMS;'+ProteinName+';'+AAname+';'+str(swapNum)+';'+str(pos) ,'Amino acid name': AAname, 'Swap size': swapNum, 'Swap start position': pos, 'Variant protein': variant, 'Variant DNA': variantDNA, 'Category': 'DMS mutant'})

# Make double mutants
pairs = ''
if args.pairwise==True:
	pairs = '_Pairwise'
	print('Generating pairwise scan of mutants')
	for AAname1, AAswap1 in AAdict.items():
		for AAname2, AAswap2 in AAdict.items():
			for pos1 in range(0, len(WTseq_AA)-1):
				for pos2 in range(pos1+1, len(WTseq_AA)):
					if pos1==pos2:
						continue
					seq = list(WTseq_AA) # set to original
					if seq[pos1]==AAswap1 or seq[pos2]==AAswap2:
						continue
					seq[pos1] = AAswap1
					seq[pos2] = AAswap2
					variant = ''.join(seq) # protein sequence
					if variant==WTseq_AA:
						continue
					dna_sequence = reverse_translate(variant)
					variantDNA = optimizeOligo(dna_sequence, BsmBIpattern)
					rows_list.append({'label':'Pairwise;'+ProteinName+';'+AAname1+';'+str(pos1)+';'+AAname2+';'+str(pos2) ,'Amino acid name': AAname1+';'+AAname2, 'Swap size': 2, 'Swap start position': str(pos1)+';'+str(pos2), 'Variant protein': variant, 'Variant DNA': variantDNA, 'Category': 'Pairwise mutant'})


# Make df
DMSdf = pd.DataFrame(rows_list, columns = ['label', 'Amino acid name','Swap size','Swap start position','Variant protein','Variant DNA', 'Category'])

print('Saving DMS variants. # Oligos:', DMSdf.shape[0])
print(str(len(set(DMSdf['Variant protein']))) + ' unique protein variants')
print(str(len(set(DMSdf['Variant DNA']))) + ' unique oligo sequences (check: should be all of them!)')

DMSdf.to_csv(args.name + 'DMS_swapUpTo' + str(args.swaps) + pairs + '.csv', index = False)

