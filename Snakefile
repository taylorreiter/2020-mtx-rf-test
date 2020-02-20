##########################
## IMPORTANT NOTE
##########################

# This snakefile is built on intermediate
# outputs from github.com/taylorreiter/ibd
# The rules that generate these intermediate
# outputs are not included in this snakefile,
# and include preprocessing and sourmash
# signature generation.
# This Snakefile is a proof-of-concept for
# correlating k-mers with quantitative/metabolite
# measurements.

import pandas as pd
import feather
from sourmash import signature
import glob
import os
from collections import Counter

m = pd.read_csv("inputs/hmp2_metadata_metatranscriptomics.tsv", sep = "\t", header = 0)
SAMPLES = m['External.ID'].unique().tolist()

rule all:
    input:
        "outputs/hash_tables/normalized_abund_hashes_wide.feather",

rule filter_signatures_to_vita_hashes_mtx:
    input:
        filt_sig = "outputs/vita_rf/vita_vars.sig",
        sigs = "outputs/sigs/{mtx}.scaled2k.sig"
    output: "outputs/filt_sigs_vita_mtx/{mtx}_filt_vita.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_vita_filtered_sigs:
    input: "outputs/filt_sigs_vita_mtx/{mtx}_filt_vita.sig"
    output: "outputs/filt_sigs_vita_named_mtx/{mtx}_filt_vita_named.sig"
    conda: "sourmash.yml"
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.mtx}_filt_vita
    '''
    
rule convert_greater_than_1_signatures_to_csv_mtx:
    input: "outputs/filt_sigs_vita_named_mtx/{mtx}_filt_vita_named.sig"
    output: "outputs/filt_sigs_vita_named_csv_mtx/{mtx}_filt_vita_named.csv"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input:
        expand("outputs/filt_sigs_vita_named_csv_mtx/{mtx}_filt_vita_named.csv", mtx = SAMPLES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/normalized_hash_abund_long.R"

rule make_hash_abund_table_wide:
    input: "outputs/hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    run:
        import pandas as pd
        import feather

        ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
        ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
        ibd_wide = ibd_wide.fillna(0)
        ibd_wide['sample'] = ibd_wide.index
        ibd_wide = ibd_wide.reset_index(drop=True)
        ibd_wide.columns = ibd_wide.columns.astype(str)
        ibd_wide.to_feather(str(output))
