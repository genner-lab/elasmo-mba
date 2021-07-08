# elasmo-mba
English Channel elasmobranch diversity paper

```bash
### SET UP PROJ ###
# download pipeline
git clone https://github.com/genner-lab/meta-fish-pipe.git
cd meta-fish-pipe
git checkout v1.0
Rscript -e "renv::restore()"
```

```r
### RETRIEVE REFERENCE LIBRARIES ###

# grab the latest reference library (run in R)
library("here")
library("tidyverse")
library("ape")
# load remote references and scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")
# read in special seqs
locals <- read_csv(file="assets/local-12s.csv")
# write out combined
reflib.orig %>% bind_rows(locals) %>% write_csv(file=here("meta-fish-lib/assets/meta-fish-lib-v243.csv"))
```

```
# add refseq
git clone https://github.com/genner-lab/refseq-reflib.git
cd refseq-reflib
mkdir temp references
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
scripts/download.sh
scripts/extract.R -p elas02
scripts/annotate.R -s 42 -p elas02



# add all the reflibs etc
# assets/contaminants-exclude-may2021.csv
# assets/meta-fish-lib-v243.csv
cp ../refseq-reflib/references/refseq206-annotated-elas02.csv assets/refseq206-annotated-elas02.csv

# just once
scripts/session-info.sh  -r assets/refseq206-annotated-elas02.csv -c assets/meta-fish-lib-v243.csv


### FOR BOTH LIBRARIES ###

# get date and time
then=$(date)

# prep
scripts/prepare-libraries.sh -p elas02 -l lib1
scripts/prepare-libraries.sh -p elas02 -l lib2
scripts/prepare-libraries.sh -p elas02 -l lib3

# make symlinks LIB1
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib1_R1.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib1/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib1_R2.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib1/fastq/R2.fastq.gz
# make symlinks LIB2
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib2_R1.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib2/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib2_R2.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib2/fastq/R2.fastq.gz
# make symlinks LIB3
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib3_R1.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib3/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/SeaDNA_Elas02_Lib3_R2.fastq.gz ~/Projects/genner-lab/elasmo-mba/meta-fish-pipe/temp/processing/elas02-lib3/fastq/R2.fastq.gz

# generate barcodes ELAS02
#F GTTGGTHAATCTCGTGCCAGC [21]
#R CATAGTAGGGTATCTAATCCTAGTTTG [27]
scripts/generate-barcodes.R -p elas02 -l lib1 -f 21 -r 27 -m assets/sequencing-master-elas02.csv
scripts/generate-barcodes.R -p elas02 -l lib2 -f 21 -r 27 -m assets/sequencing-master-elas02.csv
scripts/generate-barcodes.R -p elas02 -l lib3 -f 21 -r 27 -m assets/sequencing-master-elas02.csv


# demux
scripts/demultiplex.sh -p elas02 -l lib1 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21
scripts/demultiplex.sh -p elas02 -l lib2 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21
scripts/demultiplex.sh -p elas02 -l lib3 -f GTTGGTHAATCTCGTGCCAGC -r CATAGTAGGGTATCTAATCCTAGTTTG -t 8 -m 21


# denoise with dada2
scripts/dada2.R -p elas02 -l lib1
scripts/dada2.R -p elas02 -l lib2
scripts/dada2.R -p elas02 -l lib3

# generate stats
scripts/generate-stats.sh -p elas02 -l lib1 -t 8
scripts/generate-stats.sh -p elas02 -l lib2 -t 8
scripts/generate-stats.sh -p elas02 -l lib3 -t 8


############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############
############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############


# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p elas02

# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude-may2021.csv





```