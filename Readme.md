# Workflow Description

1.  Downloaded two source files, 
    1. One from effectorDB: https://github.com/stewart-lab/Heliobacter/blob/master/Data/uniprot\_hecliobacter\_ids.csv . This has a row for all 1500+ genes, so I chose to filter using the "T4SS effector prediction (T4SEpre)" column, selecting the 83 genes which have the value "YES."
    2. One from uniprot, https://github.com/stewart-lab/Heliobacter/blob/master/Data/uniprot-helicobacter%2Bpylori-filtered-organism\_\_Helicobacter%2Bpylori%2B(strain--.tab This is used to get the uniprot gene id, after joining the two tables on the HP\_xxxx id.
2. Looked up AA sequence with uniprot web api, e.g.  https://www.uniprot.org/uniprot/P12345.fasta
3. Translated AA to NA using the BioPython package (Bio.Data.CodonTable.generic\_by\_name["Bacterial"].back\_table ) Perhaps its important to note that this mapping is not accurate because of the many to one mapping from nucleotides to AA. (But because we edit these in the next step anyway, I thought perhaps its acceptable.)
4. Step through each protein for tiled sequences, length of 240 bp with a step size of 30 bp. 
5. Use DNA chisel to solve an optimization problem with the following parameters: 

```
def optimizeOligo(
    dna_sequence, pattern=dc.EnzymeSitePattern("BsmBI"), species="h_sapiens"
):
    """ Repurposd from Josh Tycko 2018 Script """
    problem = dc.DnaOptimizationProblem(
        sequence=dna_sequence,
        constraints=[
            dc.AvoidPattern(pattern),
            dc.EnforceGCContent(mini=0.20, maxi=0.75, window=50),
            dc.EnforceTranslation(),
        ],
        objectives=[dc.CodonOptimize(species=species)],
    )
```
