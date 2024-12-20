
#---Packages Used---------------------------------------------------------------

library(msa) #BiocManager::install("msa")
library(ape)
library(seqinr)
library(Biostrings)
library(scales)
library(lubridate)

library(ggplot2)
library(ggfortify)
library(ggtree) #BiocManager::install("ggtree")
library(rmarkdown)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(knitr)
library(reshape2)
library(tidyverse)
library(data.table)

library(pheatmap)
library(ggdendro)

library(readxl)
library(textshape)
library(growthcurver)
library(viridis)
library(gt)
library(zoo) 
library(broom)
    library(modelr)
    library(rstatix)
library(nls.multstart)


#---Global options for Markdown-------------------------------------------------

options(digits = 3)

knitr::opts_chunk$set(
    comment = "#>",
    collapse = TRUE,
    echo = FALSE,
    message = FALSE,
    warning = FALSE,
    cache = FALSE,
    #out.width = "80%",
    rows.print = 10,
    max.print = 1000,
    fig.align = "center",
    fig.width = 8,
    fig.asp = 0.618,
    fig.show = "hold"
)

options(DT.options = list(
    pageLength = 5,
    dom = 'itlp', 
    scrollX=TRUE, 
    scrollY=TRUE, 
    font_size=10
))

options(dplyr.print_min = 6, dplyr.print_max = 6)

#---Source code and input files-------------------------------------------------

source("src/lib/flat_violin.R")
source("src/s0_utilities_KB.R")
source("src/s2_barcodeProcessing_KB.R")
source("src/s3_compPlots.R")

source("https://raw.github.com/kevinsblake/GrowthRateR/main/growthrate.R")

#---Load data-------------------------------------------------------------------
analysis.file <- "data/TDase-DEGrADE_analyses_v3.xlsx"

# General
strains.summ <- read_excel(analysis.file, "Tab00_Summary", skip=1)
tcr.timeline <- read_excel(analysis.file, "Tab00_timeline", skip=1)

# MIC results
MIC.results <- read_excel(analysis.file, "Tab05_MICs")
GR.results <- read_excel(analysis.file, "Tab06_GRs_v2") 
GR.rep.results <- read_excel(analysis.file, "Tab07_GR-rep_v2") 
MIC.dosage.gen1 <- read_excel(analysis.file, "Tab10_dosage-gen1")
MIC.dosage.gen2 <- read_excel(analysis.file, "Tab10_dosage-gen2")
MIC.dosage.gen3 <- read_excel(analysis.file, "Tab10_dosage-gen3")
MIC.dosage <- rbind(MIC.dosage.gen1, MIC.dosage.gen2, MIC.dosage.gen3)

# Trees
EFF.seqs <- readAAStringSet("data/sequences/EFF_AAs_v1.fasta")
RPP.seqs <- readAAStringSet("data/sequences/RPP_AAs_v1.fasta")
DES_1.seqs <- readAAStringSet("data/sequences/DES-1_AAs_v1.fasta")
DES_2.seqs <- readAAStringSet("data/sequences/DES-2_AAs_v1.fasta")

EFF.metadata <- read_excel(analysis.file, "Tab01_EFF") # na="NA"
RPP.metadata <- read_excel(analysis.file, "Tab02_RPP")
DES_1.metadata <- read_excel(analysis.file, "Tab03_DES-1")
DES_2.metadata <- read_excel(analysis.file, "Tab04_DES-2")

# Growth rates
EXP20220602.raw <- read_excel("data/plateReader/EXP20220602/EXP20220602_v2.xlsx")
    EXP20220602.platemap <- read_excel("data/plateReader/EXP20220602/EXP20220602_platemap.xlsx")
EXP20220605.raw <- read_excel("data/plateReader/EXP20220605/EXP20220605_v2.xlsx")
    EXP20220605.platemap <- read_excel("data/plateReader/EXP20220605/EXP20220605_platemap.xlsx")
EXP20220626.raw <- read_excel("data/plateReader/EXP20220626/EXP20220626_v2.xlsx")
    EXP20220626.platemap <- read_excel("data/plateReader/EXP20220626/EXP20220626_platemap.xlsx")

    
### Directory structure --------------------------------------------------------

cur.date <- format(Sys.Date(), "%y%m%d")
fln.prfx <- stringr::str_interp("${cur.date}")
main.dir <- file.path(getwd(), "reports")

## Provide the dir name to create under main dir:
out.path <- file.path(main.dir, fln.prfx)

if (!dir.exists(out.path)){
    dir.create(out.path)
} else {
    print("output dir already exists!")
}

