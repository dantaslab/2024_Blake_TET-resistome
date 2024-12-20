
### MIC heatmaps ---------------------------------------------------------------
# remove unnecessary columns
MIC.results.subs <- subset(MIC.results, select=-c(code, mech, TET_MIC, OXY_MIC, CTC_MIC, DOX_MIC, MIN_MIC, MET_MIC, TIG_MIC, OMA_MIC, ERA_MIC))

col_order <- c("gene", "CTC", "OXY", "TET", "MIN", "MET", "DOX", "TIG", "OMA", "ERA")

MIC.results.subs <- MIC.results.subs[, col_order]
MIC.results.subs2 <- column_to_rownames(MIC.results.subs, "gene")
MIC.results.mat <- as.matrix(MIC.results.subs2)

# Create MECHANISM annotation df
annot.row <- data.frame("mech" = MIC.results$mech)
rownames(annot.row) <- MIC.results$gene

# # GEN 1 & 2 & 3
# ALL.heatmap.fname <- stringr::str_interp("${out.path}/${fln.prfx}-genes-pheatmap-ALL_v4.pdf")    
# p.MIC.ALL <- pheatmap(MIC.results.mat,
#                       color=colorRampPalette(c("#01353D", "#58B3C7","#B8FCFC"))(9),
#                       cluster_col=FALSE,
#                       annotation_row = annot.row,
#                       #cluster_rows = mat_cluster_rows,
#                       filename=ALL.heatmap.fname,
#                       annotation_colors=mech.pal2,
#                       cellheight=7.5, cellwidth=7.5,
#                       fontsize=6) 

# # GEN 1 & 2
# MIC.results.mat.1.2 <- subset(MIC.results.mat, select=-c(TIG, ERA, OMA))
# ONE.TWO.heatmap.fname <- stringr::str_interp("${out.path}/${fln.prfx}-genes-pheatmap-ONE-TWO_v1.pdf")
# 
# p.MIC.1.2 <- pheatmap(MIC.results.mat.1.2,
#                       color=colorRampPalette(c("#01353D", "#58B3C7", "#B8FCFC"))(9),
#                       #cluster_rows=FALSE,
#                       cluster_col=FALSE,
#                       #annotation_col = annot_col.df,
#                       annotation_row = annot.row,
#                       #cluster_rows = mat_cluster_rows,
#                       annotation_colors=mech.pal2,
#                       #legend=FALSE,
#                       filename=ONE.TWO.heatmap.fname,
#                       cellheight=12, cellwidth=12) +
#                    theme_pub()
# 
# # GEN 1
# MIC.results.mat.1 <- subset(MIC.results.mat.1.2, select=-c(DOX, MIN, MET))
# ONE.heatmap.fname <- stringr::str_interp("${out.path}/${fln.prfx}-genes-pheatmap-ONE_v1.pdf")
# 
# p.MIC.1.2 <- pheatmap(MIC.results.mat.1,
#                       cluster_col=FALSE,
#                       color=colorRampPalette(c("#01353D", "#58B3C7", "#B8FCFC"))(9),
#                       #annotation_col = annot_col.df,
#                       annotation_row = annot.row,
#                       #cluster_rows = mat_cluster_rows,
#                       annotation_colors=mech.pal2,
#                       filename=ONE.heatmap.fname,
#                       cellheight=12, cellwidth=12) +
#                    theme_pub()


### MICs PCA -------------------------------------------------------------------

MIC.summ <- MIC.results.subs %>%
  inner_join(strains.summ, by = c("gene" = "gene"))

# All drugs
x <- MIC.summ[2:10]
pc <- prcomp(x, scale. = TRUE) 

df <- cbind(pc$x[,1:2], MIC.summ[,11]) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(MIC.summ)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(MIC.summ)))
df$mech <- as.factor(df$mech)

p.MICs.PCA.all <- ggplot(df, aes(PC1, PC2, colour = mech)) +
  theme_pub() +
  theme(legend.position = "none",
        ) +
  geom_point() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$mech == "DES1" | df$mech == "DES2" | df$mech == "EFF" | df$mech == "RPP",]) +
  scale_color_manual(values=mech.pal) +
  xlab("PC1 (59.4%)") + ylab("PC2 (19.7%)")

# Just gen 1
y <- MIC.summ[2:4]
pc <- prcomp(y, scale. = TRUE) 

df <- cbind(pc$x[,1:2], MIC.summ[,11]) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(MIC.summ)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(MIC.summ)))
df$mech <- as.factor(df$mech)

p.MICs.PCA.gen1 <- ggplot(df, aes(PC1, PC2, colour = mech)) +
  theme_pub() +
  theme(legend.position = "none",
  ) +
  geom_point() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$mech == "DES1" | df$mech == "DES2" | df$mech == "EFF" | df$mech == "RPP",]) +
  scale_color_manual(values=mech.pal) +
  xlab("PC1 (59.4%)") + ylab("PC2 (19.7%)")

# Gen 1-2
z <- MIC.summ[2:7]
pc <- prcomp(z, scale. = TRUE) 

df <- cbind(pc$x[,1:2], MIC.summ[,11]) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(MIC.summ)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(MIC.summ)))
df$mech <- as.factor(df$mech)

p.MICs.PCA.gen2 <- ggplot(df, aes(PC1, PC2, colour = mech)) +
  theme_pub() +
  theme(legend.position = "none",
  ) +
  geom_point() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$mech == "DES1" | df$mech == "DES2" | df$mech == "EFF" | df$mech == "RPP",]) +
  scale_color_manual(values=mech.pal) +
  xlab("PC1 (70.2%)") + ylab("PC2 (19.9%)")

p.MICs.PCA.merge <- ggarrange(p.MICs.PCA.gen1, p.MICs.PCA.gen2, p.MICs.PCA.all, nrow=1, ncol=3)

p.MICs.PCA.merge

# ### Autoplot method
# MIC.pca.plot <- autoplot(pc,
#                          data = MIC.summ,
#                          colour = "mech")
# 
# MIC.pca.plot



### Growth Rates ---------------------------------------------------------------

strains.summ$code <- as.character(strains.summ$code)
GR.results$code <- as.character(GR.results$code)
GR.rep.results$code <- as.character(GR.rep.results$code)

GR.summ <- GR.results %>%
  inner_join(strains.summ, by = c("code"="code", "gene"="gene"))
GR.summ$code <- as.character(GR.summ$code)

GR.rep.summ <- GR.rep.results %>%
  inner_join(strains.summ, by = c("code"="code", "gene"="gene"))  
GR.rep.summ$code <- as.character(GR.rep.summ$code)


# GR of ea strain
p.strain.gr <- ggplot(GR.summ, aes(x=reorder(gene, -Average), y=Average, color=mech)) +
  geom_errorbar(data=GR.summ, mapping=aes(ymin=Average-StDev, ymax=Average+StDev), width=0.2, size=0.5) +
  theme_plot() +
  theme(
    legend.position = "none" 
    #aspect.ratio=1
  ) +
  ylim(0.6, 0.9) +
  scale_color_manual(values=mech.pal) +
  geom_point() +
  xlab("strain") + ylab("growth rate (h-1)")


# Bin by mechanism
## Mech bin stats
# stat.mech <- aov(Average ~ mech, GR.summ) %>% tukey_hsd()
# stat.mech.subs <- subset(stat.mech, group1 == "none" | group2 == "none")
# stat.mech.subs2 <- subset(stat.mech.subs, p.adj.signif == "*")

## GRs binned by mechanism
p.mech.gr <- ggplot(GR.summ, aes(x=reorder(mech, -Average), y=Average)) +
  theme_plot() +
  theme(
    legend.position = "none" 
    #aspect.ratio=1
  ) +
  geom_boxplot(aes(fill=mech), outlier.shape=NA) +
  scale_fill_manual(values=mech.pal) +
  #stat_pvalue_manual(stat.mech.subs2, label="p.adj.signif", y.position=0.85) +
  #stat_pvalue_manual(stat.mech.subs, label="p.adj.signif", y.position=c(0.92, 0.88, 0.86, 0.90)) +
  #geom_point() +
  geom_jitter(width=0.1, alpha=0.4) +
  ylim(0.6, 0.9) +
  xlab("mechanism") + ylab("growth rate (h-1)")

# mech.bin.fname <- stringr::str_interp("${out.path}/${fln.prfx}-mech-bin-stats_v2.pdf")
# ggsave(plot.TEST.gr.4, filename = mech.bin.fname, device=cairo_pdf, width=2.25, height=1.75)


# Plot ALL technical replicates in BIN
## Mech bin stats
# stat.rep.mech <- aov(slope ~ mech, GR.rep.summ) %>% tukey_hsd()
# stat.rep.mech.subs <- subset(stat.rep.mech, group1 == "none" | group2 == "none")
# stat.rep.mech.subs2 <- subset(stat.rep.mech.subs, p.adj.signif == "*" | p.adj.signif == "**")

# ## Plot
# plot.TEST.gr.5 <- ggplot(GR.rep.summ, aes(x=reorder(mech, -slope), y=slope)) +
#                 #geom_errorbar(data=comp.MICs.gr, mapping=aes(ymin=Average-StDev, ymax=Average+StDev), width=0.2, size=0.5) +
#                 theme_plot() +
#                 theme(legend.position = "none", 
#                     aspect.ratio=1,
#                     axis.line.x=element_blank(),
#                     axis.line.y=element_blank(),
#                     panel.border = element_rect(colour = "black", fill=NA, size=1),
#                     axis.ticks = element_line(color = "black", size=0.75),
#                     axis.text.x = element_text(angle=90, hjust=0.95, vjust=0.2)
#                     ) +
#                 #ylim(0.6, 0.9) +
#                 geom_boxplot(aes(fill=mech), outlier.shape=NA) +
#                 scale_fill_manual(values=mech.pal) +
#                 #stat_pvalue_manual(stat.rep.mech.subs2, label="p.adj.signif", y.position=c(0.92, 0.90)) +
#                 #stat_pvalue_manual(stat.rep.mech.subs, label="p.adj.signif", y.position=c(0.96, 0.92, 0.90, 0.94)) +
#                 geom_jitter(width=0.12) +
#                 xlab("mechanism") + ylab("growth rate")
# 
# print(plot.TEST.gr.5)
# 
# mech.bin.fname <- stringr::str_interp("${out.path}/${fln.prfx}-mech-bin-rep-stats_v2.pdf")
# ggsave(plot.TEST.gr.5, filename = mech.bin.fname, device=cairo_pdf, width=2) #width=5 


### MERGE GR + MICs ------------------------------------------------------------

gene.order <- c("tetX2", "tetX3", "tetX6", "tetX7", "tetX8", "tetX12", 
                "tetM", "tetO", "tetS", "tetW", "tet32", "tet36",
                "tetA", "tetB", "tetE", "tetG", "tetL", "tet39",
                "tet47", "tet50", "tet51", "tet54", "tet55", "tet56",
                "empty")

# remove unnecessary columns
comp.MICs <- subset(MIC.results, select=-c(TET, OXY, CTC, DOX, MIN, MET, TIG, OMA, ERA))
col_order <- c("gene", "code", "mech", "CTC_MIC", "OXY_MIC", "TET_MIC", "MIN_MIC", "MET_MIC", "DOX_MIC", "TIG_MIC", "OMA_MIC", "ERA_MIC")

comp.MICs <- comp.MICs[, col_order]

#melt
comp.MICs.melt <- reshape2::melt(comp.MICs, id=c("code", "mech", "gene"))
comp.MICs.melt$code <- as.character(comp.MICs.melt$code)

GR.results$code <- as.character(GR.results$code)
GR.rep.results$code <- as.character(GR.rep.results$code)

# Merge with growth rate data
comp.MICs.gr <- comp.MICs.melt %>%
  inner_join(GR.results, by="code")
comp.MICs.gr$value <- as.factor(comp.MICs.gr$value)

comp.MICs.rep.gr <- comp.MICs.melt %>%
  inner_join(GR.rep.results, by="code")
comp.MICs.rep.gr$value <- as.factor(comp.MICs.rep.gr$value)

comp.MICs.gr$gene.x <- factor(comp.MICs.gr$gene.x, levels = gene.order)

# BY MECHANISM
p.MICs.gr <- ggplot(comp.MICs.gr, aes(x=value, y=Average, color=mech)) + #gene.x
  #geom_errorbar(data=comp.MICs.gr, mapping=aes(ymin=Average-StDev, ymax=Average+StDev), width=0.2, size=0.5) +
  theme_plot() +
  theme(legend.position = "none",
        aspect.ratio = 1
  ) +
  #facet_wrap(~ variable + mech) +
  facet_wrap(~ variable) +
  scale_x_discrete(limits = c("0.25", "0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256", "512")) +
  ylim(0.6, 0.9) +
  scale_color_manual(values=mech.pal) +
  xlab ("MIC (μg/mL)") + ylab ("growth rate (h-1)") +
  geom_point()

# tmp.fname.MIC.GR <- stringr::str_interp("${out.path}/${fln.prfx}-MIC.GR_v3_KSB.pdf")
# ggsave(plot.MICs.gr, filename = tmp.fname.MIC.GR, device=cairo_pdf, width=3, height=4)

# # BY GENE
# plot.MICs.gene.gr <- ggplot(comp.MICs.gr, aes(x=value, y=Average, color=gene.x)) + #gene.x
#                 #geom_errorbar(data=comp.MICs.gr, mapping=aes(ymin=Average-StDev, ymax=Average+StDev), width=0.2, size=0.5) +
#                 theme_pub() +
#                 theme(legend.position = "right",
#                     aspect.ratio = 1,
#                     panel.border = element_rect(colour = "black", fill=NA, size=1),
#                     axis.ticks = element_line(color = "black", size=0.75),
#                     axis.text.x   = element_text(angle=90, hjust=0.95, vjust=0.2),
#                     ) +
#                 #facet_wrap(~ variable + mech) +
#                 facet_wrap(~ variable) +
#                 scale_x_discrete(limits = c("0.25", "0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256", "512")) +
#                 ylim(0.6, 0.9) +
#                 scale_color_manual(values=gene.pal) +
#                 xlab ("MIC (μg/mL)") + ylab ("growth rate (h-1)") +
#                 geom_point()
# 
# print(plot.MICs.gene.gr)
# 
# tmp.fname.MIC.gene.GR <- stringr::str_interp("${out.path}/${fln.prfx}-MIC.GR-gene_v2_KSB.pdf")
# ggsave(plot.MICs.gene.gr, filename = tmp.fname.MIC.gene.GR, device=cairo_pdf, width = 10)




### MIC dosage curves ----------------------------------------------------------

# Remove rows where abx darkened media, artificially raising OD
MIC.dosage.fil <- MIC.dosage %>%
  filter(!(abx=="TET" & conc > 256)) %>%
  filter(!(abx=="CTC" & conc > 256)) %>%
  filter(!(abx=="OXY" & conc > 256)) %>%
  filter(!(abx=="DOX" & conc > 256)) %>%
  filter(!(abx=="MET" & conc > 128)) %>%
  filter(!(abx=="MIN" & conc > 256)) %>%
  #filter(!(abx=="TIG" & conc > 256)) 
  filter(!(abx=="OMA" & conc > 32)) 
  #filter(!(abx=="ERA" & conc > 256)) 

# Wrangle
MIC.dosage.fil$conc <- as.factor(MIC.dosage.fil$conc)

# Plot
p.MIC.dose <- ggplot(MIC.dosage.fil, aes(x=conc, y=OD600, group=gene, color=gene)) +
  theme_plot() +
  theme(legend.position = "none",
        aspect.ratio = 0.75
  ) +
  facet_wrap (~ mech + abx, ncol=9) +
  geom_point(size=1) +
  geom_smooth(size=0.75,
    #method=drm,
    se=F
  ) +
  scale_color_manual(values = gene.pal) +
  geom_hline(yintercept=0.1, color="gray70") +
  xlab("drug concentration (ug/mL)") + ylab("OD600")

print(p.MIC.dose)

