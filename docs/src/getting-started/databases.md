# Databases

Some databases must be downloaded by you; the others are downloaded by the pipeline when no path is given and can be saved in `<outdir>/Databases/` for the next run.

## Overview

| Database      | Analysis Type      | Size (approx.) | Auto-setup Supported | Manual Setup/Links | Notes |
|---------------|-------------------|----------------|----------------------|--------------------|-------|
| **hostile**   | Read-based        | Variable       | Yes (default: human-t2t-hla) | [human genome refs](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/) | Bowtie2 index required |
| **sourmash**  | Read-based        | Variable      | No                   | [ctbrown’s farm](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/) | RocksDB format recommended |
| **YACHT**     | Read-based        | Variable      | No                   | [Zenodo](https://zenodo.org/communities/yacht/records), [Direct link](https://oucruaap-my.sharepoint.com/:f:/g/personal/quantnm_oucru_org/Eu87eLI5JThKg4fNcJEYXAkBUMOSsJB3_MHtNukqVs4oag?e=4tlG0b) | Edit JSON manifest after download |
| **CARD**      | Read-based/ARG detection     | ~14 GB         | Yes (default: CARD + CARD-R database)                  | [CARD download](https://card.mcmaster.ca/download) | Use `rgi_bwt` for setup |
| **geNomad**   | Assembly-based    | ~2 GB          | Yes                  | [geNomad docs](https://github.com/apcamargo/genomad#downloading-the-database) | Use `genomad download-database .` |
| **skani**     | Assembly-based    | ~50 GB         | No                   | [skani wiki](https://github.com/bluenote-1577/skani/wiki/Pre%E2%80%90sketched-databases) | Build with `skani sketch` |
| **BUSCO**     | Assembly-based    | ~60 GB         | No                   | [BUSCO docs](https://busco.ezlab.org/busco_userguide.html#download-and-automated-update) | Large output, use only if needed |
| **CheckM**    | Assembly-based    | ~1.5 GB        | Yes                  | [CheckM download](https://data.ace.uq.edu.au/public/CheckM_databases/) | Use `--save_checkm_data` to save database |
| **CheckM2**   | Assembly-based    | ~3 GB          | Yes                  | [CheckM2 download](https://zenodo.org/records/14897628) | Use `--save_checkm2_data` to save database |
| **GUNC**      | Assembly-based    | ~12.6 GB       | Yes                  | [GUNC docs](https://grp-bork.embl-community.io/gunc/usage.html#gunc-download-db) | Use `--gunc_database_type` to select |
| **Bakta**     | Assembly-based/Annotation        | 4–85 GB        | Yes (default: full)                  | [Bakta download](https://zenodo.org/records/14916843) | Use `--annotation_bakta_db_downloadtype` |
| **sylph**     | Read-based (`--readbased_taxotools sylph`) | Variable | No | [sylph pre-built databases](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases) | `--sylph_database`; sylph-tax metadata for `--sylph_taxonomy` |
| **GTDB-Tk**   | Assembly-based/Bin classification (`--bin_classifier gtdbtk`) | > 100 GB | No | [GTDB-Tk reference data](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data) | Extract the archive yourself; pass the folder with `--gtdbtk_db` |
| **AMRFinderPlus, CARD (RGI main), ResFinder/PointFinder/DisinFinder, DeepARG, VirulenceFinder, Paratype** | Bacterial genome analysis | ~5 GB (DeepARG) to ~12 GB (CARD) | Yes | See [section 13](#bacterial-genome-analysis-databases) | Saved to `${outdir}/Databases/<Tool>/` by default (`--save_<tool>_db`) |


## Setting up each database

Below are step-by-step guides for manual setup of each supported database. Follow these if you need full control or the pipeline does not support auto-download.

### hostile (host genome filtering)

- **Auto-download:** Default is `human-t2t-hla`. The pipeline can download automatically, but the database is not reusable across runs.
- **Manual build:**
  - Build a Bowtie2 index from your reference sequences:
    ```bash
    bowtie2-build path/to/sequences.fa ./<db-folder>/<index-name> --threads ${threads}
    ```
  - Download human genome references:
    ```bash
    wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/{object_name}
    ```

### sourmash

- **Option 1:** Use pre-built databases from [ctbrown’s farm](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/).
    -   [GTDB-RS226](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/)
    -   [GenBank (04/2025)](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/)
- **Option 2:** Build your own database (recommended format: RocksDB):
    ```bash
    # Create sketches from FASTA files (default: k-mer = 31, scale = 1000)
    sourmash sketch dna *.fa -o xyz.zip
    # Index sketches into a RocksDB database
    sourmash scripts index xyz.zip -o xyz.rocksdb -c ${threads} -k ${ksize}
    # Or
    sourmash index -F rocksdb xyz.rocksdb xyz.zip
    # Verify the database
    sourmash scripts check xyz.rocksdb
    ```
    > **Notes:**  
    > - Use GTDB-RS226 "genomic representatives" for best performance.  
    > - `k=31` for species-level, `k=51` for strain-level resolution.

### YACHT

- **Option 1**: Download pre-built databases from [Zenodo](https://zenodo.org/communities/yacht/records?q=&l=list&p=1&s=10&sort=newest) (GenBank 2022.03, GTDB-RS214)
    -   A direct link to the pretrained GTDB-RS226 database (ANI threshold 0.995) is available [here](https://oucruaap-my.sharepoint.com/:f:/g/personal/quantnm_oucru_org/Eu87eLI5JThKg4fNcJEYXAkBUMOSsJB3_MHtNukqVs4oag?e=4tlG0b).
- **Option 2:** Build a custom database:
    ```bash
    yacht train --ref_file xyz.zip --ksize ${ksize} --num_threads ${threads} --ani_thresh 0.995 --prefix 'xyz_ani_thresh_0.995' --outdir ./
    ```
    > **Note:** ANI threshold 0.995 recommended ([reference](https://doi.org/10.1093/bioinformatics/btae047)).
    - After download/train, edit the JSON config in the trained folder to update paths:
      ```json
      {
          "manifest_file_path": "/path/to/yacht/{database_name}/{database_name}_processed_manifest.tsv",
          "intermediate_files_dir": "/path/to/yacht/{database_name}/{database_name}_intermediate_files",
          "scale": 1000,
          "ksize": 31,
          "ani_thresh": 0.995
      }
      ```

### CARD database (for `rgi_bwt`)

- **Option 1:** Let the pipeline setup with `--enable_rgi_bwt` (do not use `--rgi_preparecarddb_dir`).
- **Option 2:** Manual setup ([RGI docs](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst#running-rgi-bwt-with-fastq-files---all-model-types)):
    ```bash
    mkdir -p ./rgi_db
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data -C ./rgi_db/
    cd ./rgi_db
    rgi card_annotation -i ./card.json > card_annotation.log 2>&1

    # Extract database version
    DB_VERSION=$(ls card_database_*_all.fasta | sed "s/card_database_v\([0-9].*[0-9]\).*/\1/")

    # Download and extract wildcard data
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
    mkdir -p wildcard
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    gunzip wildcard/*.gz
    rgi wildcard_annotation -i wildcard --card_json ./card.json -v ${DB_VERSION} > wildcard_annotation.log 2>&1

    # Finalize database (without k-mer module)
    rgi load --card_json ./card.json \
             --wildcard_annotation wildcard_database_v${DB_VERSION}.fasta \
             --wildcard_annotation_all_models wildcard_database_v${DB_VERSION}_all.fasta \
             --wildcard_index ./wildcard/index-for-model-sequences.txt \
             --wildcard_version ${DB_VERSION} \
             --card_annotation card_database_v${DB_VERSION}.fasta \
             --card_annotation_all_models card_database_v${DB_VERSION}_all.fasta \
             --debug --local 

    # Finalize database (with k-mer module)
    rgi load --card_json ./card.json \
             --wildcard_annotation wildcard_database_v${DB_VERSION}.fasta \
             --wildcard_annotation_all_models wildcard_database_v${DB_VERSION}_all.fasta \
             --wildcard_index ./wildcard/index-for-model-sequences.txt \
             --wildcard_version ${DB_VERSION} \
             --card_annotation card_database_v${DB_VERSION}.fasta \
             --card_annotation_all_models card_database_v${DB_VERSION}_all.fasta \
             --amr_kmers ./wildcard/all_amr_61mers.txt \
             --kmer_database ./wildcard/61_kmer_db.json \
             --kmer_size 61 \
             --debug --local
    ```

### geNomad

- **Option 1:** Let the pipeline setup (do not use `--genomad_db`). Use `--save_genomad_db` to save database.
- **Option 2:** Manual download ([geNomad docs](https://github.com/apcamargo/genomad#downloading-the-database)):
    ```bash
    genomad download-database .
    ```

### skani

- **Option 1:** Use pre-built databases from [skani wiki](https://github.com/bluenote-1577/skani/wiki/Pre%E2%80%90sketched-databases).
- **Option 2:** Build your own:
    ```bash
    # Option 1
    skani sketch {reference_genomes_folder}/* -o ${database_name}
    # Option 2
    find {reference_genomes_folder}/ | grep -E '\.(fna|fa|fasta)$' > ref_genomes.txt
    skani sketch -l ref_genomes.txt -o ${database_name} -t ${threads}
    ```

### BUSCO

- Download and setup as per [BUSCO documentation](https://busco.ezlab.org/busco_userguide.html#download-and-automated-update).
    > **Note:** BUSCO databases are large and outputs consume significant disk space. Use only if needed, especially for eukaryote genomes.

### CheckM/CheckM2

- **Option 1:** Let the pipeline setup (do not use `--checkm_db`/`--checkm2_db`). Use `--save_checkm_data`/`--save_checkm2_data` to save database.
- **Option 2:** Manual download:
    - [CheckM](https://data.ace.uq.edu.au/public/CheckM_databases/): ~1.5GB after extraction.
    - [CheckM2](https://zenodo.org/records/14897628): ~3GB after extraction.

### GUNC

- **Option 1:** Let the pipeline setup (do not use `--gunc_db`). Use `--gunc_database_type` to select database type (‘progenomes’ (default), ‘gtdb’) and `--gunc_save_db` to save database.
- **Option 2:** Manual download ([GUNC docs](https://grp-bork.embl-community.io/gunc/usage.html#gunc-download-db)).

### Bakta

- **Option 1:** Let the pipeline setup (do not use `--annotation_bakta_db`). Use `--annotation_bakta_db_downloadtype` to select database type (‘full’ (default), ‘light’) and `--save_bakta_db` to save database.
- **Option 2:** Manual download via [Zenodo](https://zenodo.org/records/14916843).
    - ~4GB (light) or ~85GB (full) after extraction.

### sylph and sylph-tax

- Not downloaded by the pipeline. Download a pre-built database (`*.syldb` or `*.syl2db`, e.g. GTDB R226) from the [sylph wiki](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases) and pass it with `--sylph_database` (comma-separated for several).
- For taxonomy (`--enable_sylphtax`), download the sylph-tax metadata file matching each database (e.g. `gtdb_r226_metadata.tsv.gz`, see [sylph-tax](https://github.com/bluenote-1577/sylph-tax)) and pass it with `--sylph_taxonomy`, one file per database in the same order.

### GTDB-Tk

- Not downloaded by the pipeline (the reference data is over 100 GB). Download and extract the GTDB-Tk reference data (release 226 for GTDB-Tk 2.7.2) following the [GTDB-Tk docs](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data), then pass the extracted folder with `--gtdbtk_db`. Only needed with `--bin_classifier gtdbtk`.

### Bacterial genome analysis databases

Used by `--genome_input` and `--enable_bacterial_downstream`. Every database parameter is optional:

| Tool | Parameter | When not set |
|---|---|---|
| AMRFinderPlus | `--amrfinderplus_db` | current release downloaded (`amrfinder_update`) |
| RGI main (CARD) | `--rgi_preparecarddb_dir` | CARD downloaded and prepared (shared with `rgi_bwt`) |
| ResFinder / PointFinder / DisinFinder | `--resfinder_db`, `--pointfinder_db`, `--disinfinder_db` | downloaded at `--resfinder_db_revision` |
| DeepARG | `--deeparg_db` | downloaded (`--deeparg_db_version`) |
| VirulenceFinder | `--virulencefinder_db` | downloaded at `--virulencefinder_db_revision` |
| Paratype | `--paratype_bundle` | downloaded at `--paratype_revision` |
| PlasmidFinder, ABRicate, MLST, Kleborate, SISTR, SeqSero2, Mykrobe | - | databases bundled with the containers |

A downloaded database is saved to `${outdir}/Databases/<Tool>/` (`--save_amrfinderplus_db`, `--save_rgi_db`, `--save_resfinder_db`, `--save_deeparg_db`, `--save_virulencefinder_db`, `--save_paratype_db`, all `true` by default), so later runs can pass it by path. Nothing is downloaded when no bacterial genome reaches the analysis.


> **Tip:** For most users, pipeline auto-setup is recommended where available. For large or custom databases, follow the manual instructions above. Always verify database paths and formats before running the pipeline.
