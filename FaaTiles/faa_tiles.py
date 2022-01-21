from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GENOME = "GCF_017821535.1_ASM1782153v1_genomic.fna"
GFF = "GCF_017821535.1_ASM1782153v1_genomic.gff"
PROTEIN = "GCF_017821535.1_ASM1782153v1_protein.faa"
TILED_OUT = "GCF_017821535.1_ASM1782153v1_tiled.fa"
TILE_BLOCK_SZ = 30
TILE_SZ = 240

REV_C = str.maketrans("ACGT", "TGCA")


def main():
    proteins = build_protein_metadata()
    add_na_seqs(proteins)
    make_tiled_fasta(proteins)


def make_tiled_fasta(proteins):
    with open(TILED_OUT, "w") as f:
        for id in proteins:
            for i, tile in enumerate(get_tiles(proteins[id]["na_seq"])):
                record = SeqRecord(
                    Seq(tile),
                    id=f"{id}_tile_{i:04d}",
                    name=f"{id}_tile_{i:04d}",
                    description=proteins[id]["description"],
                )
                SeqIO.write(record, f, "fasta")


def get_tiles(seq):
    beg = 0
    end = TILE_SZ
    seq_len = len(seq)
    while end < seq_len:
        yield seq[beg:end]
        beg += TILE_BLOCK_SZ
        end += TILE_BLOCK_SZ


def get_genome_seq():
    with open(GENOME, "r") as f:
        return "".join(f.read().splitlines()[1:])


def build_protein_metadata():
    proteins = {}
    i = 0
    for seq_record in SeqIO.parse(PROTEIN, "fasta"):
        i = i + 1
        proteins[seq_record.id] = {
            "description": seq_record.description,
            "aa_seq": str(seq_record.seq),
        }
    print(f" {i} total proteins read, {len(proteins)} total ids.")
    with open(GFF, "r") as f:
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


def add_na_seqs(proteins):
    genome_seq = get_genome_seq()
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
    main()
