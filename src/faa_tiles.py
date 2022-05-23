from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from pathlib import Path
import cmdlogtime
import pandas as pd
import sys
import requests
import dnachisel as dc

REV_C = str.maketrans("ACGT", "TGCA")
COMMAND_LINE_DEF_FILE = "./faa_tiles_cmdlinedef.txt"
DATA = Path("../Data")
EFFDB_FILE = DATA / "uniprot_hecliobacter_ids.csv"
ID_MAP_FILE = (
    DATA
    / "uniprot-helicobacter+pylori-filtered-organism__Helicobacter+pylori+(strain--.tab"
)
URL_BASE = "https://www.uniprot.org/uniprot"
OUTFILE = "Helio_Effector_tiled.fa"


def main(out_dir=DATA, step=30, oligoLen=240, **my_args):
    proteins = get_uniprot_seqs()

    with open(Path(out_dir) / OUTFILE, "w") as out:
        for index, row in proteins.iterrows():
            if not row["NA_Seq"]:
                continue
            cds = row["NA_Seq"][:-3]
            if len(cds) <= oligoLen:
                write_row(fp=out, chnk_num=0, id=index, row=row, seq=cds)
                continue
            for i, chunk in enumerate(
                slidingWindow(cds, winSize=oligoLen, step=step)
            ):
                write_row(
                    fp=out, chnk_num=i, id=index, row=row, seq="".join(chunk)
                )


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
    try:
        problem.resolve_constraints()
    except (TypeError, KeyError):
        print("optimization failed")
        return dna_sequence
    problem.optimize()
    optDNA = str(problem.sequence)
    return optDNA


def slidingWindow(sequence, winSize, step=1):
    """ Repurposd from Josh Tycko 2018 Script 
    
    Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""

    sequence = sequence[: len(sequence)]

    # Verify the inputs
    try:
        iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not (isinstance(winSize, int) and isinstance(step, int)):
        raise Exception("**ERROR** winSize and step must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception(
            "**ERROR** winSize must not be larger than sequence length."
        )

    # Pre-compute number of chunks to emit
    numOfChunks = int((len(sequence) - winSize) / step) + 1

    # Do the work
    for i in range(0, numOfChunks * step, step):
        yield sequence[i : i + winSize]
    # Add one more chunk, which goes to the very end
    yield sequence[-winSize:]


def write_row(fp, chnk_num, id, row, seq):
    fp.write(
        f"Heliobacter_Pylori_Effector {id} {row['HP_ids']} Tile {chnk_num+1}\n"
    )
    fp.write(f"{optimizeOligo(seq)}\n")


def get_uniprot_seqs():
    try:
        df = pd.read_csv("../Data/aa_df.tsv", sep="\t", index_col=0)
        return df
    except FileNotFoundError:
        pass

    df = pd.read_csv(EFFDB_FILE, sep="\t", index_col=0)
    df["HP_ids"] = [f'HP_{el.split(".")[1][2:]}' for el in df.index]
    df = effector_pred_filter(df)

    df_all = pd.read_csv(ID_MAP_FILE, sep="\t")
    df_all["HP_ids"] = df_all["Gene names"].apply(get_hp_id)

    df = pd.merge(
        df, df_all[["Entry", "HP_ids"]], on="HP_ids", how="inner"
    )
    df["AA_Seq"] = df.Entry.apply(get_uniprot_aa)
    df["NA_Seq"] = df.AA_Seq.apply(get_na_from_aa)
    df = df.set_index("Entry")
    df.to_csv("../Data/aa_df.tsv", sep="\t")
    return df


def effector_pred_filter(df):
    T4SE = "T4SS effector prediction (T4SEpre)"
    ELD = "Eukaryotic-like domain (EffectiveELD; PFxxx Z-Score)"
    return df[(df[T4SE] == "YES") | (df[ELD] != '-')]


def get_na_from_aa(aa):
    table = CodonTable.generic_by_name["Bacterial"].back_table
    ans = "".join([table[el] for el in aa])
    return ans.replace("U", "T")


def get_hp_id(name):
    for el in name.split():
        if el.startswith("HP_"):
            return el
    return "undefined"


def get_uniprot_aa(p):
    txt = requests.get(f"{URL_BASE}/{p}.fasta").text
    ans = "".join(txt.split("\n")[1:])
    try:
        for char in ans:
            assert char in "GAMLMFWQESKPVICYHRNDT"
        return ans
    except AssertionError:
        return ""



def make_tiled_fasta(proteins, outfile, sz, shift):
    with open(outfile, "w") as f:
        for id in proteins:
            for i, tile in enumerate(
                get_tiles(proteins[id]["na_seq"], sz, shift)
            ):
                record = SeqRecord(
                    Seq(tile),
                    id=f"{id}_tile_{i:04d}",
                    name=f"{id}_tile_{i:04d}",
                    description=proteins[id]["description"],
                )
                SeqIO.write(record, f, "fasta")


def get_tiles(seq, sz, shift):
    beg = 0
    end = sz
    seq_len = len(seq)
    while end < seq_len:
        yield seq[beg:end]
        beg += shift
        end += shift


def get_genome_seq(genome):
    with open(genome, "r") as f:
        return "".join(f.read().splitlines()[1:])


def build_protein_metadata(gff, protein_faa):
    proteins = {}
    i = 0
    for seq_record in SeqIO.parse(protein_faa, "fasta"):
        i = i + 1
        proteins[seq_record.id] = {
            "description": seq_record.description,
            "aa_seq": str(seq_record.seq),
        }
    print(f" {i} total proteins read, {len(proteins)} total ids.")
    with open(gff, "r") as f:
        for line in f:
            if not line.startswith("NZ"):
                continue
            fields = line.split("\t")
            if fields[2] == "CDS":
                (id, start, stop, is_reverse) = parse_cds(fields)
                if id not in proteins:
                    print(f"{id} not in faa file.")
                    continue
                proteins[id].update(
                    {"start": start, "stop": stop, "is_reverse": is_reverse}
                )
    return proteins


def add_na_seqs(proteins, genome):
    genome_seq = get_genome_seq(genome)
    duds = []
    big_duds = []
    for id in proteins:
        if proteins[id]["is_reverse"]:
            add_reverse_seq(id, proteins, genome_seq)
        else:
            add_forward_seq(id, proteins, genome_seq)
        proteins[id]["aa_seq_translated"] = AA_translate(
            proteins[id]["na_seq"]
        )
        if proteins[id]["aa_seq_translated"] != proteins[id]["aa_seq"]:
            if (
                proteins[id]["aa_seq_translated"][1:]
                != proteins[id]["aa_seq"][1:]
            ):
                big_duds.append(id)
            else:
                duds.append(id)
    for dud in big_duds:
        proteins.pop(dud)


def add_reverse_seq(id, proteins, genome_seq):
    proteins[id]["na_seq"] = genome_seq[
        proteins[id]["stop"] - 1 : proteins[id]["start"] - 2 : -1
    ].translate(REV_C)


def add_forward_seq(id, proteins, genome_seq):
    proteins[id]["na_seq"] = genome_seq[
        proteins[id]["start"] - 1 : proteins[id]["stop"]
    ]


def AA_translate(seq):

    table = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "_",
        "TAG": "_",
        "TGC": "C",
        "TGT": "C",
        "TGA": "_",
        "TGG": "W",
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += table[codon]
    return protein[:-1]


def parse_cds(fields):
    # ID=cds-WP_021172147.1;Parent=g
    id = fields[8].split(";")[0][7:]
    start = int(fields[3])
    stop = int(fields[4])
    is_reverse = fields[6] == "-"
    return (id, start, stop, is_reverse)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        (
            start_time_secs,
            pretty_start_time,
            my_args,
            logfile,
        ) = cmdlogtime.begin(COMMAND_LINE_DEF_FILE, sys.argv[0])
        main(my_args)
        cmdlogtime.end(logfile, start_time_secs)
