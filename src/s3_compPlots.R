### Quantify changes in mechanism / gene relative abundances -------------------
comp_analyze <- function(mix.show, drug.show, mech.show, xaxis.order){
  mix.drug <- mix.show %>%
    subset(drug == drug.show | drug == "none")
  
  if(missing(xaxis.order)){
    if (drug.show == "TET") {
      x.order <- TET.order
    }
    if (drug.show == "DOX") {
      x.order <- DOX.order
    }
    if (drug.show == "MIN") {
      x.order <- MIN.order
    }
    if (drug.show == "TIG") {
      x.order <- TIG.order
    }
    if(drug.show == "none") {
      x.order <- none.order
    }
    if(drug.show == "ATC") {
      x.order <- ATC.order
    }
  } else {
    x.order <- xaxis.order
  }
  
  mix.drug$drug_conc_timepoint <- factor(mix.drug$drug_conc_timepoint, levels=x.order)
  
  # Per change in mech rel ab
  mix.drug.change <- mix.drug %>%
    group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
    summarize(tot_mech = sum(count)) %>%
    mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
    group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
    mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
    ungroup() %>%
    group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
    summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
    ungroup()
  
  mix.drug.relab <- mix.drug %>%
    group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
    summarize(tot_mech = sum(count)) %>%
    mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
    ungroup() %>%
    group_by(drug_conc_timepoint, mech) %>%
    summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))
  
  # Per change in GENE rel ab
  mix.drug.change.gene <- mix.drug %>%
    group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
    summarize(tot_gene = sum(count)) %>%
    mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
    group_by(rep, gene) %>%
    mutate(change_per = ((gene_frac - gene_frac[1]) / gene_frac[1]) * 100) %>%
    ungroup() %>%
    group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
    summarize(mean_gene_per = mean(change_per), sd=sd(change_per)) %>%
    ungroup()
  
  # GENE rel ab
  mix.drug.relab.gene <- mix.drug %>%
    group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
    summarize(tot_gene = sum(count)) %>%
    mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
    ungroup() %>%
    group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
    summarize(mean_gene_frac = mean(gene_frac), sd=sd(gene_frac))
  
  if(missing(mech.show)){
    
    # Merge all into list
    mylist <- list()
    mylist[[1]] <- mix.drug
    mylist[[2]] <- mix.drug.change
    mylist[[3]] <- mix.drug.relab
    mylist[[4]] <- mix.drug.change.gene
    mylist[[5]] <- mix.drug.relab.gene
    
  } else {
    if(mech.show == "EFF") {
      mech.list <- EFF.list
    }
    if(mech.show == "RPP") {
      mech.list <- RPP.list
    }
    if(mech.show == "DES1") {
      mech.list <- DES1.list
    }
    if(mech.show == "DES2") {
      mech.list <- DES2.list
    }
    # Subset for just mech
    mix.drug.relab.gene.mech <- mix.drug.relab.gene %>%
      subset(gene %in% mech.list)
    
    mix.drug.change.gene.mech <- mix.drug.change.gene %>%
      subset(gene %in% mech.list)
    
    # Merge all into list
    mylist <- list()
    mylist[[1]] <- mix.drug
    mylist[[2]] <- mix.drug.change
    mylist[[3]] <- mix.drug.relab
    mylist[[4]] <- mix.drug.change.gene
    mylist[[5]] <- mix.drug.relab.gene
    mylist[[6]] <- mix.drug.relab.gene.mech
    mylist[[7]] <- mix.drug.change.gene.mech
    
  }

  return(mylist)
  
}


### Plot three-panel competition figures ---------------------------------------
comp_plot <- function(mix.show, drug.show, mech.show, xaxis.order){ 
  
  mylist <- comp_analyze(mix.show, drug.show, mech.show, xaxis.order) 
  
  p01a.mix.drug.barplot <- ggplot(mylist[[1]], aes(x=drug_conc_timepoint, y=count, fill=gene)) +
    geom_bar(stat="identity", position="fill") +
    theme_plot() +
    scale_fill_manual(values=gene.pal) +
    labs(y="rel. abundance")
  
  p01b.mix.drug.mech.per <- ggplot(mylist[[2]], aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
    geom_point() + 
    geom_line() +
    geom_errorbar(mylist[[2]], mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
    theme_plot() +
    ylim(NA, 350) + 
    scale_color_manual(values=mech.pal) +
    labs(y="% change")
  
  p01c.mix.drug.gene.mech <- ggplot(mylist[[7]], aes(x=drug_conc_timepoint, y=mean_gene_per, color=gene, group=gene)) +
    geom_point() + 
    geom_line() +
    geom_errorbar(mylist[[7]], mapping=aes(ymin=mean_gene_per-sd, ymax=mean_gene_per+sd), width=0.2, size=0.5) +
    theme_plot() +
    scale_color_manual(values=gene.pal) +
    #ylim(-100, NA) +
    labs(y="% change")
  
  p01.mix.drug.merge <- ggarrange(p01a.mix.drug.barplot, p01b.mix.drug.mech.per, p01c.mix.drug.gene.mech, ncol=3)
  
  return(p01.mix.drug.merge)  
  
}

### Plot two-panel competition figures (for mech-specific mixes) ---------------
comp_plot_mech <- function(mix.show, drug.show, mech.show, xaxis.order){ 
  
  mylist <- comp_analyze(mix.show, drug.show, mech.show, xaxis.order) 
  
  p01a.mix.drug.barplot <- ggplot(mylist[[1]], aes(x=drug_conc_timepoint, y=count, fill=gene)) +
    geom_bar(stat="identity", position="fill") +
    theme_plot() +
    scale_fill_manual(values=gene.pal) +
    labs(y="rel. abundance")
  
  p01c.mix.drug.gene.mech <- ggplot(mylist[[7]], aes(x=drug_conc_timepoint, y=mean_gene_per, color=gene, group=gene)) +
    geom_point() + 
    geom_line() +
    geom_errorbar(mylist[[7]], mapping=aes(ymin=mean_gene_per-sd, ymax=mean_gene_per+sd), width=0.2, size=0.5) +
    theme_plot() +
    scale_color_manual(values=gene.pal) +
    ylim(-100, NA) +
    labs(y="% change")
  
  p01.mix.drug.merge <- ggarrange(p01a.mix.drug.barplot, p01c.mix.drug.gene.mech, ncol=3) # 2 or 3
  
  return(p01.mix.drug.merge)  
  
}

comp_more_plot <- function(mix.show, drug.show, mech.show){
  
  mylist <- comp_analyze(mix.show, drug.show, mech.show) 
  
  # Specific comparisons
  mix.drug.relab.gene.stats <- mylist[[1]] %>%
    group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
    summarize(tot_gene = sum(count)) %>%
    mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
    ungroup() %>%
    group_by(gene) %>%
    pairwise_t_test(gene_frac ~ drug_conc_timepoint, p.adjust.method = "bonferroni")
  
  mix.drug.relab.gene.stats.sub <- mix.drug.relab.gene.stats %>%
    subset(p.adj.signif != "ns")
  
  p00c.mix.drug.gene.bar <- ggplot(mylist[[5]], aes(x=drug_conc_timepoint, y=mean_gene_frac, fill=gene)) +
    geom_bar(stat="identity") +
    geom_errorbar(mylist[[5]], mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2) +
    theme_plot() +
    #stat_pvalue_manual(mix.drug.relab.gene.stats.sub, label="p.adj.signif") + #, y.position=0.08
    ylim(0, NA) +
    facet_wrap(~ gene, scales="free") +
    scale_fill_manual(values=gene.pal) +
    labs(y="avg. abundance")
  
  return(p00c.mix.drug.gene.bar)
  
}





# # Quantify changes in mechanism / gene relative abundances
# ### Before option to not include "mech.show"
#
# comp_analyze <- function(mix.show, drug.show, mech.show){
#   mix.drug <- mix.show %>%
#     subset(drug == drug.show | drug == "none")
#   
#   if (drug.show == "TET") {
#     x.order <- TET.order
#   }
#   if (drug.show == "DOX") {
#     x.order <- DOX.order
#   }
#   if (drug.show == "MIN") {
#     x.order <- MIN.order
#   }
#   if (drug.show == "TIG") {
#     x.order <- TIG.order
#   }
#   if(drug.show == "none") {
#     x.order <- none.order
#   }
#   
#   if(mech.show == "EFF") {
#     mech.list <- EFF.list
#   }
#   if(mech.show == "RPP") {
#     mech.list <- RPP.list
#   }
#   if(mech.show == "DES1") {
#     mech.list <- DES1.list
#   }
#   if(mech.show == "DES2") {
#     mech.list <- DES2.list
#   }
#   
#   mix.drug$drug_conc_timepoint <- factor(mix.drug$drug_conc_timepoint, levels=x.order)
#   
#   # Per change in mech rel ab
#   mix.drug.change <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
#     summarize(tot_mech = sum(count)) %>%
#     mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
#     group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
#     mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
#     summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
#     ungroup()
#   
#   mix.drug.relab <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
#     summarize(tot_mech = sum(count)) %>%
#     mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, mech) %>%
#     summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))
#   
#   # Per change in GENE rel ab
#   mix.drug.change.gene <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
#     summarize(tot_gene = sum(count)) %>%
#     mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
#     group_by(rep, gene) %>%
#     mutate(change_per = ((gene_frac - gene_frac[1]) / gene_frac[1]) * 100) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
#     summarize(mean_gene_per = mean(change_per), sd=sd(change_per)) %>%
#     ungroup()
#   
#   # GENE rel ab
#   mix.drug.relab.gene <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
#     summarize(tot_gene = sum(count)) %>%
#     mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
#     ungroup() %>%
#     group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
#     summarize(mean_gene_frac = mean(gene_frac), sd=sd(gene_frac))
#   
#   # Subset for just mech
#   mix.drug.relab.gene.mech <- mix.drug.relab.gene %>%
#     subset(gene %in% mech.list)
#   
#   mix.drug.change.gene.mech <- mix.drug.change.gene %>%
#     subset(gene %in% mech.list)
#   
#   # Merge all into list
#   mylist <- list()
#   mylist[[1]] <- mix.drug
#   mylist[[2]] <- mix.drug.change
#   mylist[[3]] <- mix.drug.relab
#   mylist[[4]] <- mix.drug.change.gene
#   mylist[[5]] <- mix.drug.relab.gene
#   mylist[[6]] <- mix.drug.relab.gene.mech
#   mylist[[7]] <- mix.drug.change.gene.mech
#   
#   return(mylist)
#   
# }
# 
# 

# # Plot competition results
# comp_plot <- function(mix.show, drug.show, mech.list, x.order){
#   
#   mix.drug <- mix.show %>%
#     subset(drug == drug.show | drug == "none")
#   
#   mix.drug$drug_conc_timepoint <- factor(mix.drug$drug_conc_timepoint, levels=x.order)
#   
#   p01a.mix.drug.barplot <- ggplot(mix.drug, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
#     geom_bar(stat="identity", position="fill") +
#     theme_plot() +
#     scale_fill_manual(values=gene.pal) +
#     labs(y="rel. abundance")
#   
#   # Per change in mech rel ab
#   mix.drug.change <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
#     summarize(tot_mech = sum(count)) %>%
#     mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
#     group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
#     mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
#     summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
#     ungroup()
#   
#   mix.drug.relab <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
#     summarize(tot_mech = sum(count)) %>%
#     mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, mech) %>%
#     summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))
#   
#   p01b.mix.drug.mech.per <- ggplot(mix.drug.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
#     geom_point() + 
#     geom_line() +
#     geom_errorbar(mix.drug.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
#     theme_plot() +
#     ylim(-100, 350) +
#     scale_color_manual(values=mech.pal) +
#     labs(y="% change")
#   
#   # Per change in GENE rel ab
#   mix.drug.change.gene <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
#     summarize(tot_gene = sum(count)) %>%
#     mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
#     group_by(rep, gene) %>%
#     mutate(change_per = ((gene_frac - gene_frac[1]) / gene_frac[1]) * 100) %>%
#     ungroup() %>%
#     group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
#     summarize(mean_gene_per = mean(change_per), sd=sd(change_per)) %>%
#     ungroup()
#   
#   mix.drug.relab.gene <- mix.drug %>%
#     group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
#     summarize(tot_gene = sum(count)) %>%
#     mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
#     ungroup() %>%
#     group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
#     summarize(mean_gene_frac = mean(gene_frac), sd=sd(gene_frac))
#   
#   ### TURN mech.list into way to make EFF.list
#   mix.drug.relab.gene.mech <- mix.drug.relab.gene %>%
#     subset(gene %in% mech.list)
#   
#   # Subset for just mech
#   mix.drug.change.gene.mech <- mix.drug.change.gene %>%
#     subset(gene %in% mech.list)
#   
#   p01c.mix.drug.gene.mech <- ggplot(mix.drug.change.gene.mech, aes(x=drug_conc_timepoint, y=mean_gene_per, color=gene, group=gene)) +
#     geom_point() + 
#     geom_line() +
#     geom_errorbar(mix.drug.change.gene.mech, mapping=aes(ymin=mean_gene_per-sd, ymax=mean_gene_per+sd), width=0.2, size=0.5) +
#     theme_plot() +
#     scale_color_manual(values=gene.pal) +
#     ylim(-100, NA) +
#     labs(y="% change")
#   
#   p01.mix.drug.merge <- ggarrange(p01a.mix.drug.barplot, p01b.mix.drug.mech.per, p01c.mix.drug.gene.mech, ncol=3)
#   
#   return(p01.mix.drug.merge)
#   
# 
# }
