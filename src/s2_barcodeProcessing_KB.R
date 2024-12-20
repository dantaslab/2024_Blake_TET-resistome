
bc_process <- function(exp.name, output=c("all", "summ")){
  gene.order <- c("tetX2", "tetX3", "tetX6", "tetX7", "tetX8", "tetX12", 
                  "tetM", "tetO", "tetS", "tetW", "tet32", "tet36",
                  "tetA", "tetB", "tetE", "tetG", "tetL", "tet39",
                  "tet47", "tet50", "tet51", "tet54", "tet55", "tet56",
                  "empty")
  
  exp.dname <- stringr::str_interp("data/barseq/${exp.name}")
  exp.mapname <- stringr::str_interp("data/barseq/${exp.name}/mapping.csv")
  exp.mapping <- read_csv(exp.mapname)
  
  # Data wrangling + renaming
  filelist <- list.files(path = exp.dname, pattern = "*.txt", full.names=TRUE)
  samplenames <- list.files(path = exp.dname, pattern = ".txt")
  samplenames2 <- gsub("tabUMIs.txtcounts.txt", "", samplenames)
  
  datalist <- lapply(filelist, read.table, header=TRUE)
  colnames <- c("barcode", "count")
  datalist2 <- lapply(datalist, setNames, colnames)
  names(datalist2) <- samplenames2
  
  # Merge strain summary with barcode counts >> show barcodes with count <1
  bc.counts <- lapply(datalist2, merge, strains.summ, by="barcode", all=TRUE)
  
  # Merge all together, adding the listname as a column value
  bc.merge <- purrr::map_df(bc.counts, data.frame, .id='rxn')
  
  # Merge bc counts with sample metadata (mapping file)
  bc.merge.map <- left_join(bc.merge, exp.mapping, by="rxn") %>%
    unite(mix_rep, c(mix, rep), sep="-", remove=FALSE) %>%
    unite(mix_rep_timepoint, c(mix_rep, timepoint), sep="-", remove=FALSE) %>%
    unite(mix_timepoint, c(mix, timepoint), sep="-", remove=FALSE) %>%
    mutate(conc2 = as.character(conc)) %>%
    unite(drug_conc, c(drug, conc2), sep="-", remove=FALSE) %>%
    unite(conc_timepoint, c(conc2, timepoint), sep="-", remove=FALSE) %>%
    unite(drug_conc_timepoint, c(drug_conc, timepoint), sep="-", remove=FALSE) %>%
    unite(mix_drug_conc_timepoint, c(mix, drug_conc_timepoint), sep="-", remove=FALSE) %>%
    unite(mix_drug_conc_timepoint_rep, c(mix_drug_conc_timepoint, rep), sep="-", remove=FALSE) %>%
    unite(gene_mix_drug_conc_timepoint, c(gene, mix_drug_conc_timepoint), sep="-", remove=FALSE) %>%
    mutate(count = replace(count, count < 10, NA)) # MAKE AUTOMATIC >> COUNT # RXNS
  
  bc.merge.map$gene <- factor(bc.merge.map$gene, levels = gene.order)
  bc.merge.map$count[is.na(bc.merge.map$count)] <- 0
  
  # Threshold - Make zero any counts below the threshold (0.1%)
  bc.merge.filter <- bc.merge.map %>%
    group_by(rxn) %>%
    mutate(tot_count = sum(count))  %>%
    group_by(rxn, gene) %>%
    mutate(gene_frac = count/tot_count) %>%
    ungroup() %>%
    mutate(count=replace(count, gene_frac <= 0.001, 0)) %>%
    ungroup()
    
  # Summarize replicates
  bc.summary <- bc.merge.filter %>%
    group_by(mix, mix_timepoint, drug_conc_timepoint, mix_drug_conc_timepoint, gene, drug, conc2, mech) %>%
    summarize(mean = mean(count), sd=sd(count))
  
  if(missing(output)){
    output <- "all"
  }
  if(output == "all"){
    bc_output <- bc.merge.map
  }
  if(output == "summ"){
    bc_output <- bc.summary
  }
  
  return(bc_output)
}


rare_process <- function(exp.name){

  exp.dname <- stringr::str_interp("data/barseq/${exp.name}")
  exp.mapname <- stringr::str_interp("data/barseq/${exp.name}/mapping.csv")
  exp.mapping <- read_excel("data/TDase-DEGrADE_analyses_v3.xlsx", "Tab08d_val02_rare")
  
  # Data wrangling + renaming
  filelist <- list.files(path = exp.dname, pattern = "*.txt", full.names=TRUE)
  samplenames <- list.files(path = exp.dname, pattern = ".txt")
  samplenames2 <- gsub("tabUMIs.txtcounts.txt", "", samplenames)
  
  datalist <- lapply(filelist, read.table, header=TRUE)
  colnames <- c("barcode", "count")
  datalist2 <- lapply(datalist, setNames, colnames)
  names(datalist2) <- samplenames2
  
  # Merge strain summary with barcode counts >> show barcodes with count <1
  bc.counts <- lapply(datalist2, merge, strains.summ, by=c("barcode"), all=TRUE)
  
  # Merge all together, adding the listname as a column value
  bc.merge <- purrr::map_df(bc.counts, data.frame, .id='rxn')
  
  # Merge bc counts with sample metadata (mapping file)
  bc.merge.map <- left_join(bc.merge, exp.mapping, by=c("barcode", "rxn")) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    group_by(rxn) %>%
    mutate(measab=round(count/sum(count), 5)) %>% 
    ungroup
  
  return(bc.merge.map)
}
