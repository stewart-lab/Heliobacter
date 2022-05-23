import faa_tiles
import pandas as pd
import requests
import xmltodict

PROTEIN_DAT_FILE = faa_tiles.DATA / "UP000000429_85962.dat"
PFAM_URL_BASE = "https://pfam.xfam.org/family"


def main():
    protein_df = get_effector_proteins_with_pfam_data()
    pfam_id_set = get_pfam_set(protein_df)
    print(
        "By intersecting the effector_db list with the UP429_85962.dat file,\n",
        f"\twe find {protein_df.shape[0]} out of 112 are there. Of those, \n",
        f"\t{protein_df[protein_df['Pfam_csv']==''].shape[0]} have no annotated",
        f"Pfam domains. There are {len(pfam_id_set)} unique Pfam_ids represented.",
    )
    pfam_df = build_pfam_df(pfam_id_set)
    pfam_df.to_csv(faa_tiles.DATA / "pfam_data.csv")
    protein_df.to_csv(faa_tiles.DATA / "protein_df.csv")


def build_pfam_df(pfam_id_set):
    pfam_dict = {el: get_pfam_sub_dict(el) for el in pfam_id_set if el}
    return pd.DataFrame.from_dict(pfam_dict, orient="index")


def get_pfam_sub_dict(pfam_id):
    txt = requests.get(f"{PFAM_URL_BASE}/{pfam_id}?output=xml").text
    my_dict = xmltodict.parse(txt)["pfam"]["entry"]
    comment = ""
    if "comment" in my_dict:
        comment = my_dict["comment"]
    return {
        "Pfam_id": my_dict["@id"],
        "description": my_dict["description"],
        "comment": comment,
    }


def get_pfam_set(protein_df):
    pfam_ids = [el.split(",") for el in protein_df["Pfam_csv"]]
    return {el for row in pfam_ids for el in row}


def get_effector_proteins_with_pfam_data():
    proteins = faa_tiles.get_uniprot_seqs()
    return proteins.merge(
        get_protein_dat_df(), how="left", left_index=True, right_index=True
    )


def get_protein_dat_df():
    with open(PROTEIN_DAT_FILE, "r") as f:
        data = {k: v for (k, v) in get_prot_records(f)}
        return pd.DataFrame([data]).T.rename(columns={0: "Pfam_csv"})


def get_prot_records(f):
    # Key is in record's 2nd line, so can use first for whileloop.
    while f.readline():
        k = get_key(f.readline())
        pfams = get_pfam_list_and_goto_end_of_record(f)
        yield (k, ",".join(pfams))


def get_key(line):
    # 'AC   O05732;\n'
    assert line.split()[0] == "AC"
    assert line[-2] == ";"
    return line.strip().split()[1][:-1]


def get_pfam_list_and_goto_end_of_record(f):
    pfams = list()
    line = f.readline()
    while not line.startswith("//"):
        if line.startswith("DR   Pfam;"):
            pfams.append(line.split(";")[1].strip())
        line = f.readline()
    return pfams


if __name__ == "__main__":
    main()
