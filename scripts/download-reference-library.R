#!/usr/bin/env Rscript

# libs
suppressMessages({
    library("here")
    library("tidyverse")
    library("ape")
})

# load remote references and scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

# read in special seqs
locals <- suppressMessages(read_csv(file=here("assets/local-12s.csv")))

# write out combined
reflib.orig %>% 
    bind_rows(locals) %>% 
    write_csv(file=here("meta-fish-pipe/assets/meta-fish-lib-v243.csv"))

# report
writeLines("\nReference library saved to 'meta-fish-pipe/assets/meta-fish-lib-v243.csv'\n")
