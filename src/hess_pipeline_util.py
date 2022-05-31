import swalign
import re
import pandas as pd
from pathlib import Path
# from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import socket
import pdb

HOSTNAME = socket.gethostname()
if "chtc" in HOSTNAME:
    DATA = Path("../Data")
    FASTQ_HOME = Path("/staging/steill")
else:
    DATA = Path.home() / "Hess/Data"
    FASTQ_HOME = DATA / "Hess_Fastqs"
REF_EDIT_SITE = "gtgcaccttgaagcgcatgaactccttgatgatggccatgttatcctcctcgcccttgctcaccattgggccaggattctcctcgacatc"
DEFAULT_GLOB_PATTERN = "Nuc_*.fastq.gz"
DEFAULT_DOMAIN_FILE = FASTQ_HOME / "NucDomain_oligos.csv"  # DDRTile_oligos.csv
# https://pypi.org/project/swalign/
# choose your own values here… 2 and -1 are common.
MATCH = 2
MISMATCH = -1
SW = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(MATCH, MISMATCH))
DOMAINS = {}
N_D_MAX = 10  # number of domains in R2 (dev only)


def main_pipeline():
    df_result = pd.DataFrame()
    for r1, r2 in get_pairs(n=10):
        domain_i = get_best_domain(r2)
        cigar_j = get_cigar(r1)
        bump_count(df_result, domain_i, cigar_j)
    print(df_result.fillna(0))


def bump_count(df, i, j):
    try:
        df.loc[i, j] += 1
    except Exception:
        df.loc[i, j] = 1


def build_domain_dict(domain_file):
    return {k: seq for (k, seq) in get_keys_seqs(domain_file)}


def get_keys_seqs(domain_file):
    with open(domain_file, "r") as f:
        for line in f:
            try:
                pieces = re.split(",|;", line)
                assert len(pieces) > 4
                seq = pieces[-1].upper().strip()
                key = "_".join(pieces[:-1])
                yield (key, seq)
            except AssertionError:
                pass


def get_best_domain(r2, domain_file=DEFAULT_DOMAIN_FILE, n_d_max=N_D_MAX):
    if not DOMAINS:
        DOMAINS.update(build_domain_dict(domain_file))
    best_score = -99999
    cnt = 0
    for key, seq in DOMAINS.items():
        cnt += 1
        if not cnt % 100:
            print(f"{cnt/len(DOMAINS)*100:.2f}% Done...")
        alignment = SW.align(seq, r2)
        if alignment.score > best_score:
            best_score = alignment.score
            best_key = key
        if cnt >= n_d_max:
            break

    return best_key


def get_cigar(r1):
    alignment = SW.align(REF_EDIT_SITE, r1)
    return alignment.cigar_str

def get_paired_files():
    all_files = FASTQ_HOME.glob(DEFAULT_GLOB_PATTERN)
    common_headers = {el.name.split(".")[0][:-1] for el in all_files}
    for h in sorted(list(common_headers)):
        yield (FASTQ_HOME / f"{h}1.fastq.gz", FASTQ_HOME / f"{h}2.fastq.gz")


def get_pairs(n=10, glob_ptrn=DEFAULT_GLOB_PATTERN):
    for (r1_file, r2_file) in get_paired_files():
        count = 0
        with gzip.open(r1_file, "rt") as fh_r1, gzip.open(
            r2_file, "rt"
        ) as fh_r2:
            r1_iter = FastqGeneralIterator(fh_r1)
            r2_iter = FastqGeneralIterator(fh_r2)
            r1, r2 = r1_iter.__next__(), r2_iter.__next__()
            yield (r1[1].strip(), r2[1].strip())
            count = +1
            if count == n:
                continue


def sequence_strider(glob_ptrn, seqs_per_file=10, step_size=100):
    for ko_file in FASTQ_HOME.glob(glob_ptrn):
        cnt = 0
        with gzip.open(ko_file, "rt") as f:
            for rec in FastqGeneralIterator(f):
                cnt += 1
                if not cnt % step_size:
                    yield rec[1].lower()
                if cnt >= seqs_per_file * step_size:
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
    tgt = "CTGAACCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCATTGGGCCAGGATTCTCCTCGACAT"
    tgt = tgt.lower()
    assert get_first_indel_pos(REF_EDIT_SITE, tgt) == 62
    assert get_first_indel_pos(REF_EDIT_SITE, REF_EDIT_SITE) is None


if __name__ == "__main__":
    main_pipeline()
