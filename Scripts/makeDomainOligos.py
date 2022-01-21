'''
Make Pfam domain library

Input: database download from UniProt
	can query for things like localization, certain domains, species, etc

Output: library of oligos to express those short domains in the HT-recruit vector

October 2018
Josh Tycko
'''

from prody import *
import pandas as pd
import numpy as np
import math
from itertools import islice
import argparse
import sys
from dnachisel import *

parser = argparse.ArgumentParser(description='Get sequence of best transcript')
parser.add_argument('in_file', type=str,
                    help='Input file from UniProt site')
# parser.add_argument('domain_file', type=str,
#                     help='Input domain file of genes from the UniProt table list')
parser.add_argument('out_file', type=str,
                    help='Output oligo filename start')
parser.add_argument('residues', type=int,
                    help='Number of residues per element')

# Optional arguments:
parser.add_argument('-o', '--optimize', dest='opt',
                    help='Optimize the DNA sequences',
                    type=str, default='')
parser.add_argument('-s', dest = 'step', type=str,
                    help='sliding window step size', default = 'scan')
parser.add_argument('--swaps', dest='swaps', type=int,
                    help='Max number of consecutive AA swaps (eg. 3 gives single, double, and triple swaps', default = 0)
parser.add_argument('-c', '--controls', dest = 'controls', type=int,
                    help='Number of random control sequences to include',
                    default = 500)
args = parser.parse_args()

###############################################################################
def optimizeOligo(dna_sequence, pattern):
	problem = DnaOptimizationProblem(
			sequence=dna_sequence,
			constraints=[AvoidPattern(pattern),
			             EnforceGCContent(mini=0.20, maxi=0.75, window=50),
			             EnforceTranslation()],
			objectives=[CodonOptimize(species='h_sapiens')]
			) # should use , mode = 'harmonized'
	try:
		problem.resolve_constraints()
	except (TypeError, KeyError):
		print('optimization failed')
		return dna_sequence
	problem.optimize()
	optDNA = str(problem.sequence)
	return optDNA

###############################################################################

# Parse Database
# read in input uniProtTable
uniProtTable = pd.read_table(args.in_file)
domainFile = args.in_file[:-3] + '_domains.csv'

try:
	domainDf = pd.read_csv(domainFile)
except IOError:
	domainDf = pd.DataFrame(columns = ['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names', 'Organism', 'Protein length', 'Annotation', 'Keywords', 'Protein existence', 'Subcellular location [CC]', 'Gene ontology (biological process)', 'Sequence similarities', 'Protein families', 'Protein Sequence', 'Domain Sequence', 'evalue', 'start', 'end', 'Domain length', 'Domain accession', 'Domain ID', 'Domain type'])

gene_col = 'Entry name'

# get list of genes that needs to be filled out
gene_listA = [gene for gene in list(set(uniProtTable.loc[uniProtTable['Cross-reference (Pfam)'].notnull()][gene_col].tolist())) if str(gene) != "nan"]


# If missing Pfam domains, find their locations
gene_listB = []
for gene in gene_listA:
	if gene not in set(domainDf[gene_col].tolist()):
		gene_listB.append(gene)

print('Looking up Pfam domains for ' + str(len(gene_listB)) + ' genes out of ' + str(len(uniProtTable.index)))

gene_count = 0
for gene in gene_listB:
	gene_count += 1
	try:
		query = searchPfam(gene) # gets a dict from Pfam
	except ValueError:
		try:
			query = searchPfam(uniProtTable.loc[uniProtTable[gene_col] == gene, 'Entry'].values[0])
		except ValueError:
			print('Could not find on UniProt with gene name or UniProt ID')
			continue
	if len(query) < 1:
		continue
	rows_list = []
	for domainAcc, values in query.iteritems():
		# Get the domain information
		domType = values['type']
		domID = values['id']

		# Iterate over locations of that domain
		for loc in range(len(values['locations'])):
			start = int(values['locations'][loc]['start'])
			end = int(values['locations'][loc]['end'])
			evalue = float(values['locations'][loc]['evalue'])
			domSeq = uniProtTable.loc[uniProtTable[gene_col]==gene]['Sequence'].str[start-1: end-1].tolist()[0] #10/19 this was getting different sequence than uniprot website
				#reason is that python 0-indexes and Uniprot 1-indexes so we should look from start-1 to end-1
				#This actually made a difference for MPP8
				#Fixed it for next time design libraries with this code...
				#Good news is we were almost always extending the domain

			rows_list.append({'Domain Sequence': domSeq, 'evalue': evalue, 'start': start, 'end': end, 'Domain length': end-start, 'Domain accession': domainAcc, 'Domain ID': domID, 'Domain type': domType})

	# Copy the protein information
	genedf = pd.DataFrame(rows_list, columns = ['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names', 'Organism', 'Protein length', 'Annotation', 'Keywords', 'Protein existence', 'Subcellular location [CC]', 'Gene ontology (biological process)', 'Sequence similarities', 'Protein families', 'Protein Sequence', 'Domain Sequence', 'evalue', 'start', 'end', 'Domain length', 'Domain accession', 'Domain ID', 'Domain type'])

	genedf.loc[:,['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names', 'Organism', 'Protein length', 'Annotation', 'Keywords', 'Protein existence', 'Subcellular location [CC]', 'Gene ontology (biological process)', 'Sequence similarities', 'Protein families', 'Protein Sequence']] = uniProtTable.loc[uniProtTable[gene_col] == gene][['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names', 'Organism', 'Length', 'Annotation', 'Keywords', 'Protein existence', 'Subcellular location [CC]', 'Gene ontology (biological process)', 'Sequence similarities', 'Protein families', 'Sequence']].values.tolist()[0]

	domainDf = domainDf.append(genedf, ignore_index = True)
	if gene_count % 100 == 0:
		print('Saving domain file' + domainFile)
		domainDf.to_csv(domainFile, index = False)

# Save the final domain file
print('Saving domain file' + domainFile)
domainDf.to_csv(domainFile, index = False)

# Optional: also include "Regions" - these are not Pfam domains but they are annotated with a function
# note - we should have done this!

# Find the amino acid chain we want to express
rows_list = []
noDomainFound = 0
ProteinTooShort = 0
DomainTooLong = 0
for index, row in domainDf.iterrows():
	start = int(row['start'])
	end = int(row['end'])
	domSeq = row['Domain Sequence']
	ProtSeq = row['Protein Sequence']
	proteinID = row['Entry name']
	domID = row['Domain ID']
	domLen = len(domSeq)

	# Filter out domains we don't want
	if domID in ['zf-C2H2']: # about half of the short domains, by number
		continue

	# Confirm we are matching the domain with the right full-length protein
	if domSeq not in ProtSeq:
		noDomainFound += 1
		continue

	if len(ProtSeq) < args.residues:
		ProteinTooShort += 1
		continue

	if domLen > args.residues:
		DomainTooLong += 1
		continue

	if domLen <= args.residues:
		if start - (args.residues - domLen)/2 < 2:
			start = 1 # don't use Start M
			end = args.residues + 1
		elif end + (args.residues - domLen)/2 > len(ProtSeq):
			end = len(ProtSeq)
			start = end - args.residues
		else:
			start = start - (args.residues - domLen)/2
			end = end + (args.residues - domLen)/2
	if end - start < args.residues:
		start -= 1

	ExtDomSeq = ProtSeq[start:end]
	
	# Find DNA sequence using dnachisel
	ExtDNASeq = reverse_translate(ExtDomSeq)

	# Codon optimize it, adjust GC, and remove BsmBI sites
	if args.opt == 'CodonOpt':
		ExtDNASeq = optimizeOligo(ExtDNASeq, enzyme_pattern("BsmBI"))
		
	rows_list.append({'Gene entry name': proteinID, 'Domain ID': domID, 'Domain start': start, 'Domain length': domLen, 'Protein sequence': ProtSeq, 'Domain sequence': domSeq, 'Extended Domain sequence': ExtDomSeq, 'Extended Domain DNA sequence': ExtDNASeq})

oligodf = pd.DataFrame(rows_list, columns = ['Gene entry name', 'Domain ID', 'Domain start', 'Domain length', 'Protein sequence', 'Domain sequence', 'Extended Domain sequence', 'Extended Domain DNA sequence'])

# Optional: Alanine scan across it


# Remove redundant domains
dups = sum(oligodf['Domain sequence'].duplicated())
print('Duplicates', dups)
oligodf['Gene entry names of redundant domains'] = None
for domSeq in set(list(oligodf[oligodf['Domain sequence'].duplicated(keep=False)]['Domain sequence'])): # this makes it way faster
			proteinIDs = list(oligodf.loc[oligodf['Domain sequence'] == domSeq]['Gene entry name'])
			if len(set(proteinIDs)) > 1:
				oligodf.loc[oligodf['Domain sequence'] == domSeq, 'Gene entry names of redundant domains'] = str(proteinIDs)
oligodf = oligodf.drop_duplicates('Domain sequence')

# Handle few bug entries in oligodf from Titin (no oligo was found)


# Save oligos
print(oligodf.shape[0])
oligodf.to_csv(args.out_file + 'less' + str(args.residues) + 'aa_' + args.opt + 'withoutZFC2H2.csv', index = False)





