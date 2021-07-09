# Environmental DNA captures elasmobranch diversity in a temperate marine ecosystem

Code and data for:

Liu, Z., Collins, R.A., Baillie, C., Rainbird, S., Griffiths, A.M., Sims, D.W., Mariani, S. & Genner, M.J. (2021). Environmental DNA captures elasmobranch diversity in a temperate marine ecosystem. _Insert Journal_. [https://doi.org/xxx](https://doi.org/xxx).

---

1. Clone this project repository onto your system and install all required R package versions.

```bash
# download repo and R packages
git clone https://github.com/genner-lab/elasmo-mba.git
cd elasmo-mba
Rscript -e "renv::restore()"
```

2. Download the meta-fish-pipe (v1.0) bioinformatics module. Be sure to read instructions and install all system software as instructed at [github.com/genner-lab/meta-fish-pipe](https://github.com/genner-lab/meta-fish-pipe).

```bash
# download and install pipeline
git clone https://github.com/genner-lab/meta-fish-pipe.git
cd meta-fish-pipe
git checkout v1.0
Rscript -e "renv::restore()"
cd ..
```

3. Download the custom UK fish reference library (v243) from [github.com/genner-lab/meta-fish-lib](https://github.com/genner-lab/meta-fish-lib).

```bash
# download the custom uk fish reference library
scripts/download-reference-library.R
```

4. Download the NCBI RefSeq reference library (v206) from [github.com/genner-lab/refseq-reflib](https://github.com/genner-lab/refseq-reflib).

```bash
# download RefSeq reference library
git clone https://github.com/genner-lab/refseq-reflib.git
cd refseq-reflib
git checkout v1.0
Rscript -e "renv::restore()"
mkdir temp references
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
scripts/download.sh
scripts/extract.R -p elas02
scripts/annotate.R -s 42 -p elas02
rm temp/duckdb
cd ..
```

5. Copy sample files into the meta-fish-pipe module.

```bash
# copy across files to bioinformatics pipeline
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv
cp assets/sequencing-master-elas02.csv meta-fish-pipe/assets/sequencing-master-elas02.csv
cp refseq-reflib/references/refseq206-annotated-elas02.csv meta-fish-pipe/assets/refseq206-annotated-elas02.csv
```

6. Generate session information for the bioinformatics pipeline.

```bash
# generate session stats
cd meta-fish-pipe
scripts/session-info.sh -r assets/refseq206-annotated-elas02.csv -c assets/meta-fish-lib-v243.csv
```

7. Prepare the directories for each library.

```bash 
# prepare each library
scripts/prepare-libraries.sh -p elas02 -l lib1
scripts/prepare-libraries.sh -p elas02 -l lib2
scripts/prepare-libraries.sh -p elas02 -l lib3
```

8. Copy across symbolic links to the fastq files.

```bash
# make symlinks LIB1
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib1_R1.fastq.gz temp/processing/elas02-lib1/fastq/R1.fastq.gz
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib1_R2.fastq.gz temp/processing/elas02-lib1/fastq/R2.fastq.gz
# make symlinks LIB2
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib2_R1.fastq.gz temp/processing/elas02-lib2/fastq/R1.fastq.gz
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib2_R2.fastq.gz temp/processing/elas02-lib2/fastq/R2.fastq.gz
# make symlinks LIB3
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib3_R1.fastq.gz temp/processing/elas02-lib3/fastq/R1.fastq.gz
ln -s -r /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib3_R2.fastq.gz temp/processing/elas02-lib3/fastq/R2.fastq.gz
```

9. Generate barcode files for each library.

```bash
# generate barcodes files for ELAS02
scripts/generate-barcodes.R -p elas02 -l lib1 -f 21 -r 27 -m assets/sequencing-master-elas02.csv
scripts/generate-barcodes.R -p elas02 -l lib2 -f 21 -r 27 -m assets/sequencing-master-elas02.csv
scripts/generate-barcodes.R -p elas02 -l lib3 -f 21 -r 27 -m assets/sequencing-master-elas02.csv
```

10. Demultiplex each library with cutadapt.

```bash
# demultiplex with cutadapt
scripts/demultiplex.sh -p elas02 -l lib1 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21
scripts/demultiplex.sh -p elas02 -l lib2 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21
scripts/demultiplex.sh -p elas02 -l lib3 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21
```

11. Denoise each library with dada2.

```bash
# denoise with dada2
scripts/dada2.R -p elas02 -l lib1
scripts/dada2.R -p elas02 -l lib2
scripts/dada2.R -p elas02 -l lib3
```

12. Generate pipeline output stats.

```bash
# generate bioinfomatics stats
scripts/generate-stats.sh -p elas02 -l lib1 -t 8
scripts/generate-stats.sh -p elas02 -l lib2 -t 8
scripts/generate-stats.sh -p elas02 -l lib3 -t 8
```

13. Run the taxonomic assignement steps and combine results from each library 

```bash
# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p elas02
# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude.csv
cd ..
```

14. Run the statistical analyses (this step should preferably be run line-by-line in the preferred R console to assess results in turn).  

```bash
# run stats
scripts/stats-figures.R
```
