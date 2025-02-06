import csv
import pandas as pd
import numpy as np
import urllib
import gzip
import screed

basename = "vmr_MSL39_v4"
out_dir = f"output.{basename}"
logs_dir = os.path.join(out_dir, 'logs')

configfile: "conf/db-config.yml"
ASSEMBLY_FASTA_DIR = "genbank/genomes"
CURATED_FASTA_DIR = f"genbank/curated/{basename}" # curated fasta files

# full VMR
# first use find-assembly-accessions.py to generate this file.
vmr_file = 'inputs/VMR_MSL39.v4_20241106.acc.tsv'
vmr = pd.read_csv(vmr_file, sep='\t')

sourmash_params = config['sourmash_params']
moltypes = sourmash_params.keys()

# expand files for rule all:
param_combos = []
for moltype in moltypes:
    ksizes = sourmash_params[moltype]["ksize"]
    scaleds = sourmash_params[moltype]["scaled"]
    combo = expand(f"{moltype}-k{{k}}-sc{{sc}}", k=ksizes, sc=scaleds)
    param_combos.extend(combo)

# select as needed
#param_combos = [x for x in param_combos if not x.endswith('-sc1')]
param_combos = [x for x in param_combos if x.endswith('-sc50')]
print(param_combos)

# rules for building params for sourmash sketching
"""
build single parameter string for sourmash sketching
"""
def build_param_str(moltype, override_scaled=None):
    ksizes = sourmash_params[moltype]['ksize']
    if override_scaled:
        scaled = override_scaled
    else:
        scaled = min(sourmash_params[moltype]['scaled'])
    k_params = ",".join([f"k={k}" for k in ksizes])
    param_str = f"-p {moltype},{k_params},scaled={scaled},abund"
    return param_str

"""
build multiple params for all sourmash sketching
"""
def build_params(sourmash_params, override_scaled=None):
    param_str = []
    for moltype in sourmash_params.keys():
        param_str.append(build_param_str(moltype, override_scaled))
    return " ".join(param_str)


wildcard_constraints:
    acc = '[^/]+',
    vmr_acc = '[^/]+',

rule all:
    input:
        expand(os.path.join(out_dir, f"{basename}.sc{{scaled}}.zip"), scaled=[1000]), #10
        expand(os.path.join(out_dir, f"{basename}.{{params}}.rocksdb/CURRENT"), params=param_combos),
        os.path.join(out_dir, "blastn", f"{basename}.index.nhr"),


### Rules for ICTV GenBank Assemblies:
rule download_assembly_summary:
    output:
        summary = 'genbank/assembly_summary.viral.txt',
        historical = 'genbank/assembly_summary_historical.txt', # includes suppressed accessions, etc
    shell: 
        """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt > {output.summary}
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary_historical.txt > {output.historical} 
        """

rule acc_to_directsketch:
    input:
        vmr_file = vmr_file,
        good_acc = 'genbank/assembly_summary.viral.txt',
        bad_acc = 'genbank/assembly_summary_historical.txt',
    output:
        ds_csv = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
        curated_ds = os.path.join(out_dir, f"{basename}.urlsketch.csv"),
        curate_info = os.path.join(out_dir, f"{basename}.ncfasta-to-curate.csv"),
        suppressed = os.path.join(out_dir, f"{basename}.suppressed.csv"),
        lengths = os.path.join(out_dir, f"{basename}.lengths.csv"),
        lineages = os.path.join(out_dir, f"{basename}.lineages.csv"),
    shell:
        """
        python acc-to-directsketch.py \
            --vmr_file {input.vmr_file} \
            --good_acc {input.good_acc} \
            --bad_acc {input.bad_acc} \
            --ds_csv {output.ds_csv} \
            --curated_ds {output.curated_ds} \
            --curate_info {output.curate_info} \
            --suppressed {output.suppressed} \
            --lengths {output.lengths} \
            --lineages {output.lineages} \
            --basename {basename}
        """

rule directsketch_assembly_datasets:
    input:
        csvfile = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.sc1.gbsketch.zip"),
        failed = os.path.join(out_dir, f"{basename}.gbsketch-assemblies-failed.csv"),
        ch_failed = os.path.join(out_dir, f"{basename}.gbsketch-assemblies-checksum-failed.csv"),
    params:
        fastadir= ASSEMBLY_FASTA_DIR,
        param_str = build_params(sourmash_params),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "directsketch", f"{basename}.gbsketch.log")
    benchmark: os.path.join(logs_dir, "directsketch", f"{basename}.gbsketch.benchmark")
    shell:
        """
        sourmash scripts gbsketch -o {output.zipf} {input.csvfile} -n 1 \
                                  --keep-fasta --fastas {params.fastadir} \
                                  --genomes-only {params.param_str} --retry-times 5 \
                                  --failed {output.failed} --checksum-fail {output.ch_failed} 2> {log}
        """

rule download_assembly_datasets:
    input:
        csvfile = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
    output:
        failed = os.path.join(out_dir, f"{basename}.gbsketch-download-assemblies-failed.csv"),
        ch_failed = os.path.join(out_dir, f"{basename}.gbsketch-download-assemblies-checksum-failed.csv"),
    params:
        fastadir= ASSEMBLY_FASTA_DIR,
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "directsketch", f"{basename}.download-gbsketch.log")
    benchmark: os.path.join(logs_dir, "directsketch", f"{basename}.download-gbsketch.benchmark")
    conda: "directsketch-188.yml"
    shell:
        """
        sourmash scripts gbsketch {input.csvfile} -n 10 \
                                  --keep-fasta --fastas {params.fastadir} \
                                  --download-only --genomes-only --retry-times 15 \
                                  --failed {output.failed} --checksum-fail {output.ch_failed} 2> {log}
        """

# # Where we don't have assembly datasets, curate fasta files from GenBank nuccore
# here, we use directsketch to merge the component fastas into a single sketch AND take ranges where needed
rule directsketch_curated:
    input: 
        os.path.join(out_dir, f"{basename}.urlsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.sc1.urlsketch.zip"), 
        failed = os.path.join(out_dir, f"{basename}.urlsketch-failed.csv"),
    params:
        fastadir= CURATED_FASTA_DIR,
        param_str = build_params(sourmash_params),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "urlsketch", f"{basename}.urlsketch.log")
    benchmark: os.path.join(logs_dir, "urlsketch", f"{basename}.urlsketch.benchmark")
    shell:
        """
        sourmash scripts urlsketch {input} -o {output.zipf} -n 1 \
                                  --keep-fasta --fastas {params.fastadir} \
                                  --genomes-only {params.param_str} --retry-times 5 \
                                  --failed {output.failed} 2> {log}
        """

rule directsketch_curated_download:
    input: 
        os.path.join(out_dir, f"{basename}.urlsketch.csv"),
    output:
        failed = os.path.join(out_dir, f"{basename}.urlsketch-download-failed.csv"),
    params:
        fastadir= CURATED_FASTA_DIR,
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "urlsketch", f"{basename}.download.urlsketch.log")
    benchmark: os.path.join(logs_dir, "urlsketch", f"{basename}.download.urlsketch.benchmark")
    conda: "directsketch-main.yml"
    shell:
        """
        sourmash scripts urlsketch {input} -n 9 --download-only \
                                  --keep-fasta --fastas {params.fastadir} \
                                  --genomes-only --retry-times 15 \
                                  --failed {output.failed} 2> {log}
        """

rule combine_sigs:
    input:
        directsketch = os.path.join(out_dir, f"{basename}.sc1.gbsketch.zip"),
        curated = os.path.join(out_dir, f"{basename}.sc1.urlsketch.zip"),
    output:
        combined = os.path.join(out_dir, f"{basename}.sc1.zip"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "sigcat", f"{basename}.log")
    benchmark: os.path.join(logs_dir, "sigcat", f"{basename}.benchmark")
    shell:
        """
        sourmash sig cat {input.directsketch} {input.curated} -o {output.combined} 2> {log}
        """


rule downsample_sigs:
    input:
         combined = os.path.join(out_dir, f"{basename}.sc1.zip"),
    output:
         downsampled = os.path.join(out_dir, f"{basename}.sc{{scaled}}.zip"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    log: os.path.join(logs_dir, "downsample", f"{basename}.sc{{scaled}}.log")
    benchmark: os.path.join(logs_dir, "downsample", f"{basename}.sc{{scaled}}.benchmark")
    shell:
        """
        sourmash sig downsample {input.combined} -o {output} --scaled {wildcards.scaled} 2> {log}
        """

rule index_rocksdb:
    input:
        os.path.join(out_dir, f"{basename}.sc{{scaled}}.zip"),
    output:
        rocksdb_current = os.path.join(out_dir, f"{basename}.{{moltype}}-k{{ksize}}-sc{{scaled}}.rocksdb/CURRENT"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    params:
        rocksdb_dir = directory(os.path.join(out_dir, f"{basename}.{{moltype}}-k{{ksize}}-sc{{scaled}}.rocksdb")),
    log: os.path.join(logs_dir, "index", f"{basename}.{{moltype}}-k{{ksize}}-sc{{scaled}}.log")
    benchmark: os.path.join(logs_dir, "index", f"{basename}.{{moltype}}-k{{ksize}}-sc{{scaled}}.benchmark")
    shell:
        """
        sourmash scripts index -o {params.rocksdb_dir} {input} \
                         --scaled {wildcards.scaled} --moltype {wildcards.moltype} \
                         --ksize {wildcards.ksize} 2> {log}
        """


## rules for BLAST db generation

# get lengths of sequences in FASTA files + combine all into single FASTA file
rule combine_fasta_and_save_length_info:
    input:
        #zipf = os.path.join(out_dir, f"{basename}.sc1.urlsketch.zip"), # use zip to make sure we have the fasta files
        gb_fail = os.path.join(out_dir, f"{basename}.gbsketch-download-assemblies-failed.csv"),
        gb_c_failed = os.path.join(out_dir, f"{basename}.urlsketch-download-failed.csv"),
        urlsketch_csv = os.path.join(out_dir, f"{basename}.urlsketch.csv"),
        gbsketch_csv = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
    output: 
        combined= os.path.join(out_dir, f"{basename}.fna.gz"),
        lengths= os.path.join(out_dir, f"{basename}.lengths.csv"),
    params:
        curated_fasta_dir = CURATED_FASTA_DIR,
        assembly_fasta_dir = ASSEMBLY_FASTA_DIR,
    log: os.path.join(logs_dir, "combine-fasta", f"{basename}.log")
    benchmark: os.path.join(logs_dir, "combine-fasta", f"{basename}.benchmark")
    shell:
        """
        python combine-fasta.py --urlsketch-csv {input.urlsketch_csv} \
                                --gbsketch-csv {input.gbsketch_csv} \
                                --combined {output.combined} \
                                --lengths {output.lengths} \
                                --curated-fasta-dir {params.curated_fasta_dir} \
                                --assembly-fasta-dir {params.assembly_fasta_dir} 2> {log}
        """

# # Rule to build BLAST index for the combined gzipped fasta file
rule build_blast_nucl_index:
    input:
        fasta = os.path.join(out_dir, f"{basename}.fna.gz"),
    output:
        index = os.path.join(out_dir, "blastn", f"{basename}.index.nhr")
    params:
        title = os.path.join(f"{basename}"),
        out_base = os.path.join(out_dir, "blast", f"{basename}.blastn.index")
    log:  os.path.join(logs_dir, "blast", f"{basename}.blastn-index.log")
    benchmark:  os.path.join(logs_dir, "blast", f"{basename}.blastn-index.benchmark")
    conda: "conf/env/blast.yml"
    shell:
        """
        gunzip -c {input.fasta} | makeblastdb -in - -dbtype nucl -parse_seqids \
               -out {params.out_base} -title {params.title} 2> {log}
        """

# rule build_prot_index:
#     input:
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fasta = os.path.join(out_dir, "{basename}.protein.fa.gz"),
#     output:
#         index = os.path.join(out_dir, "diamond", "{basename}.protein.fa.gz" + ".dmnd"),
#     log:  os.path.join(logs_dir, "diamond-index", "{basename}.log")
#     benchmark:  os.path.join(logs_dir, "diamond-index", "{basename}.benchmark")
#     conda: "conf/env/diamond.yml"
#     shell:
#         """
#         diamond makedb --in {input.fasta} --db {output.index} --quiet 2> {log}
#         """
