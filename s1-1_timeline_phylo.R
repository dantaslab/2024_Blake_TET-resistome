
### Timeline --------------------------------------------------------------------

# Plot
timeline_plot <- ggplot(tcr.timeline,aes(x=year,y=0, col=gen, label=milestone)) +
  labs(col="Milestones") +
  scale_color_manual(values=mech.pal) +
  theme_pub() +
  theme(
    axis.line.y=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    #axis.text.x =element_blank(),
    axis.ticks.x =element_blank(),
    axis.line.x =element_blank(),
    legend.position = "none"
  ) +
  geom_hline(yintercept=0, color="black", size=0.3) +
  geom_segment(aes(y=position, yend=0,xend=year), color="black", size=0.2) +
  geom_point(aes(y=0), size=3) +
  xlim(1940, 2024) +
  geom_text(aes(y=text_position, label=milestone), size=2.5)


### EFF tree --------------------------------------------------------------------

# Generate MSA
EFF.seqs.aln <- msa(EFF.seqs, "ClustalOmega") # can also use "Muscle" etc

# Generate percentage pairwise identity matrix with ape
EFF.seqs.aln.2 <- msaConvert(EFF.seqs.aln, type="ape::AAbin") # Convert to ape-friendly format
EFF.seqs.dist <- dist.gene(EFF.seqs.aln.2, "percentage") # Create pairwise distance matrix
EFF.seqs.dist.mat <- as.matrix(EFF.seqs.dist)
EFF.seqs.dist.per <- 100 * (1 - EFF.seqs.dist.mat) # Make percentage pairwise identity

# # Export
# mat1.fname <- stringr::str_interp("${out.path}/${fln.prfx}-EFF_perID_wide_v1.csv")
# write.csv(EFF.seqs.dist.per, file=mat1.fname)

# p1 - Make base tree
EFF.tree <- nj(EFF.seqs.dist) # neighbor-joining method
p.EFF.tree <- ggtree(EFF.tree) + geom_tiplab() + xlim(0, 1)

# p2 - Base tree + colored names designating which used in competition
p2.EFF.tree <- p.EFF.tree %<+% EFF.metadata +
  geom_tiplab(aes(fill=factor(compete)),
              color = "black", geom = "label", label.padding = unit(0.15, "lines"), label.size = 0) +
  scale_fill_manual(values=c ("#FFFFFF","#EB8F94","#00BFC4")) +
  labs(title="EFF sequences") +
  theme(legend.position ="none")


### RPP tree -------------------------------------------------------------------

# Generate MSA
RPP.seqs.aln <- msa(RPP.seqs, "ClustalOmega") # can also use "Muscle" etc

# Generate percentage pairwise identity matrix with ape
RPP.seqs.aln.2 <- msaConvert(RPP.seqs.aln, type="ape::AAbin") # Convert to ape-friendly format

RPP.seqs.dist <- dist.gene(RPP.seqs.aln.2, "percentage") # Create pairwise distance matrix
RPP.seqs.dist.mat <- as.matrix(RPP.seqs.dist)
RPP.seqs.dist.per <- 100 * (1 - RPP.seqs.dist.mat) # Make percentage pairwise identity

# # Export
# mat1.fname <- stringr::str_interp("${out.path}/${fln.prfx}-RPP_perID_wide_v1.csv")
# write.csv(RPP.seqs.dist.per, file=mat1.fname)

# p1 - Make base tree
RPP.tree <- nj(RPP.seqs.dist) # neighbor-joining method
p.RPP.tree <- ggtree(RPP.tree) + geom_tiplab() + xlim(0, 1)

# p2 - Base tree + colored names designating competition
p2.RPP.tree <- p.RPP.tree %<+% RPP.metadata +
  geom_tiplab(aes(fill=factor(compete)),
              color = "black", geom = "label", label.padding = unit(0.15, "lines"), label.size = 0) +
  scale_fill_manual(values=c ("#FFFFFF","#EB8F94","#00BFC4")) +
  labs(title="RPP sequences") +
  theme(legend.position ="none")


### DES1 tree ------------------------------------------------------------------

# Generate MSA
DES_1.seqs.aln <- msa(DES_1.seqs, "ClustalOmega") # can also use "Muscle" etc

# Generate percentage pairwise identity matrix with ape
DES_1.seqs.aln.2 <- msaConvert(DES_1.seqs.aln, type="ape::AAbin") # Convert to ape-friendly format
DES_1.seqs.dist <- dist.gene(DES_1.seqs.aln.2, "percentage") # Create pairwise distance matrix
DES_1.seqs.dist.mat <- as.matrix(DES_1.seqs.dist)
DES_1.seqs.dist.per <- 100 * (1 - DES_1.seqs.dist.mat) # Make percentage pairwise identity

# # Export
# mat1.fname <- stringr::str_interp("${out.path}/${fln.prfx}-DES-1_perID_wide_v1.csv")
# write.csv(DES_1.seqs.dist.per, file=mat1.fname)

# p1 - Make base tree
DES_1.tree <- nj(DES_1.seqs.dist) # neighbor-joining method
p.DES_1.tree <- ggtree(DES_1.tree) + geom_tiplab() + xlim(0, 0.2)

# p2 - Base tree + colored names designating ordered/functional
p2.DES_1.tree <- p.DES_1.tree %<+% DES_1.metadata +
  geom_tiplab(aes(fill=factor(compete)),
              color = "black", geom = "label", label.padding = unit(0.15, "lines"), label.size = 0) +
  scale_fill_manual(values=c ("#FFFFFF","#00BFC4")) + 
  labs(title="DES1 sequences") +
  theme(legend.position ="none")


### DES2 tree ------------------------------------------------------------------

# Generate MSA
DES_2.seqs.aln <- msa(DES_2.seqs, "ClustalOmega") # can also use "Muscle" etc

# Generate percentage pairwise identity matrix with ape
DES_2.seqs.aln.2 <- msaConvert(DES_2.seqs.aln, type="ape::AAbin") # Convert to ape-friendly format
DES_2.seqs.dist <- dist.gene(DES_2.seqs.aln.2, "percentage") # Create pairwise distance matrix
DES_2.seqs.dist.mat <- as.matrix(DES_2.seqs.dist)
DES_2.seqs.dist.per <- 100 * (1 - DES_2.seqs.dist.mat) # Make percentage pairwise identity

# # Export
# mat1.fname <- stringr::str_interp("${out.path}/${fln.prfx}-DES-2_perID_wide_v1.csv")
# write.csv(DES_2.seqs.dist.per, file=mat1.fname)

# p1 - Make base tree
DES_2.tree <- nj(DES_2.seqs.dist) # neighbor-joining method
p.DES_2.tree <- ggtree(DES_2.tree) + geom_tiplab() + xlim(0, 1)

# p2 - Base tree + colored names designating ordered/functional
p2.DES_2.tree <- p.DES_2.tree %<+% DES_2.metadata +
  geom_tiplab(aes(fill=factor(compete)),
              color = "black", geom = "label", label.padding = unit(0.15, "lines"), label.size = 0) +
  scale_fill_manual(values=c ("#FFFFFF","#00BFC4")) + 
  labs(title="DES2 sequences") +
  theme(legend.position ="none")

