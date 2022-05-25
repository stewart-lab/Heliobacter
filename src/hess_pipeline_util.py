import swalign
from pathlib import Path
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import socket

HOSTNAME = socket.gethostname()
if 'chtc' in HOSTNAME:
    DATA = Path('../Data')
    FASTQ_HOME = Path("/staging/steill")
else:
    DATA = Path.home() / 'Hess/Data'
    FASTQ_HOME = DATA / 'Hess_Fastqs'
REF = 'gtgcaccttgaagcgcatgaactccttgatgatggccatgttatcctcctcgcccttgctcaccattgggccaggattctcctcgacatc'
DEFAULT_GLOB_PATTERN = "Nuc_*.fastq.gz"
DEFAULT_DOMAIN_FILE = FASTQ_HOME / "NucDomain_oligos.csv" # DDRTile_oligos.csv
# https://pypi.org/project/swalign/
# choose your own values here… 2 and -1 are common.
MATCH = 2
MISMATCH = -1
SW = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(MATCH, MISMATCH))

def main():
    pass

def get_best_domain(r2, domain_file = DEFAULT_DOMAIN_FILE):
    return str(10*int(r2))

def get_cigar(r1):
    return f'M{r1}'

def get_pairs(n=10, glob_ptrn=DEFAULT_GLOB_PATTERN):
    
    for i in range(n):
        yield (i, -i)

def sequence_strider(glob_ptrn, seqs_per_file=10, step_size=100):
    for ko_file in FASTQ_HOME.glob(glob_ptrn):
        cnt = 0
        with gzip.open(ko_file,'rt') as f:
            for rec in FastqGeneralIterator(f):
                cnt += 1
                if not cnt%step_size:
                    yield rec[1].lower()
                if cnt >= seqs_per_file*step_size: 
                    break

def get_first_indel_pos(ref, target):
    alignment = SW.align(ref, target)
    r_pos = alignment.r_pos
    for cig in alignment.cigar:
        if cig[1] not in "ID":
            r_pos += cig[0]
        else:
            return r_pos
    return None

def test_first_indel_pos():
    tgt = 'CTGAACCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCATTGGGCCAGGATTCTCCTCGACAT'
    tgt = tgt.lower()
    assert get_first_indel_pos(REF, tgt) == 62
    assert get_first_indel_pos(REF, REF) is None
    

if __name__ == "__main__":
    main()