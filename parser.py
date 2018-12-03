import os
import csv
from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number
from biothings.utils.dataload import merge_duplicate_rows, dict_sweep
from utils.hgvs import get_hgvs_from_vcf
from itertools import groupby

VALID_COLUMN_NO = 22

'''this parser is for BioMuta v3.0(BioMuta3 Complete Dataset) downloaded from
https://hive.biochemistry.gwu.edu/cgi-bin/prd/biomuta/servlet.cgi'''


# convert one snp to json
def _map_line_to_json(df):
    # specific variable treatment
    chrom = df["chr_id"]
    pos = df["chr_pos"]
    if chrom == 'M':
        chrom = 'MT'

    ref = df["ref_nt"]
    alt = df["alt_nt"]

    HGVS = get_hgvs_from_vcf(chrom, int(pos), ref, alt, mutant_type=False)
    
    transcript_id = clean_data(df["transcript_id"], ("-",))
    peptide_id = clean_data(df["peptide_id"], ("-",))
    uniprot_ac = clean_data(df["uniprot_ac"], ("-",))
    refseq_ac = clean_data(df["refseq_ac"], ("-",))
    cds_pos = clean_data(df["cds_pos"], ("-",))
    pep_pos = clean_data(df["pep_pos"], ("-",))
    uniprot_pos = clean_data(df["uniprot_pos"], ("-",))
    ref_aa = clean_data(df["ref_aa"], ("-",))
    alt_aa = clean_data(df["alt_aa"], ("-",))
    mut_freq = clean_data(df["mut_freq"], ("-",))
    data_src = clean_data(df["data_src"], ("-",))
    do_id = clean_data(df["do_id"], ("-",))
    do_name_id, do_name = do_name_split(df["do_name"])
    assert do_id == do_name_id, "do_id mismatch!"

    uberon_id = to_list(df["uberon_id"])
    gene_name = clean_data(df["gene_name"], ("-",))
    pmid_list = to_list(df["pmid_list"])
    site_prd = clean_data(df["site_prd"], ("-",))
    site_ann = df["site_ann"]


# load as json data
    one_snp_json = {
        "_id": HGVS,
        "biomuta": {
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'transcript_id': transcript_id,
            'peptide_id': peptide_id,
            'uniprot_ac': uniprot_ac,
            'refseq_ac': refseq_ac,
            'cds_pos': cds_pos,
            'pep_pos': pep_pos,
            'uniprot_pos': uniprot_pos,
            'ref_aa': ref_aa,
            'alt_aa': alt_aa,
            'mut_freq': mut_freq,
            'data_src': data_src,
            'do_id': {
                        do_id : do_id,
                        do_name : do_name
                        },
            'uberon_id': uberon_id,
            'gene_name': gene_name,
            'pmid': pmid_list,
            'site_prd': site_prd,
            'site_ann': site_ann,
        }
    }
    one_snp_json = value_convert_to_number(one_snp_json)
    return one_snp_json


def clean_index(s):
    return s.lower().replace("/", "_").replace("-", "_").replace("(", "_").replace(")", "").replace("#", "")


def clean_data(d, vals):
    if d in vals:
        return None
    else:
        return d

def to_list(s, sep=";"):
    if s:
        return s.split(sep)
    else:
        return None

def do_name_split(do_name):
    if not do_name:
        return None, None
    do_id, do_name = do_name.split(" / ")
    do_id = do_id.split(":")[1]
    return do_id, do_name


# open file, parse, pass to json mapper
def load_data(data_folder):
    input_fn = os.path.join(data_folder,"biomuta-master.csv")
    open_file = open(input_fn)
    db_biomuta = csv.reader(open_file)
    index = next(db_biomuta)
    assert len(index) == VALID_COLUMN_NO, "Expecting %s columns, but got %s" % (VALID_COLUMN_NO, len(index))
    index = [clean_index(s) for s in index]
    biomuta = (dict(zip(index, row)) for row in db_biomuta)
    biomuta = filter(lambda row: row["index"] != "", biomuta)
    json_rows = map(_map_line_to_json, biomuta)
    json_rows = (row for row in json_rows if row)
    json_rows = sorted(json_rows, key=lambda row: row["_id"])
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "biomuta") for rg in row_groups)
    return (unlist(dict_sweep(row, vals=[None, ])) for row in json_rows)
