# 2025-ictv-vmr-database

Build sourmash databases for the new ICTV release.

First, clone this repo and cd in. Then:

1. Install and activate conda env

```
mamba env create -f environment.yml
mamba activate 2025-ictv-vmr-database
```

2. If rebuilding this database, unzip the provided accession file: `gunzip inputs/VMR_MSL40.v1.20250307.acc.tsv.gz` so it will be available for the snakemake workflow. Proceed to #3.

Alternatively, if building from a new ICTV VMR file, you'll need to get the Assembly Accessions first like so
```
python scripts/get-assembly-acc2.py -i VMR_MSL40.v1.20250307.xlsx -o VMR_MSL40.v1.20250307.acc.tsv \
                                    -s "VMR MSL40" --output-directsketch
```
> If this fails to get all on the first attempt, move the output file to a tmp file, e.g. `mv VMR_MSL40.v1.20250307.acc.tsv VMR_MSL40.v1.20250307.acc1.tsv`. Then provide that existing file to the script using `--existing-vmr VMR_MSL40.v1.20250307.acc1.tsv`. You can do this any number of times, always moving the latest output file to ensure you have all existing results. 

3. Modify the config file (`conf/db-config.yml`) to point to the final `*acc.tsv` file and change signature params as desired.

4. Run the snakemake workflow:
```
snakemake -c 30 -j 1 -n 1
```
* Remove `-n` (dry run) to actually run the workflow.*
