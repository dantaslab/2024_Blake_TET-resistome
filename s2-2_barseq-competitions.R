
### Data wrangling -------------------------------------------------------------

# Wrangle t0 
EXP20230207 <- bc_process("20230207_mix0")
EXP20230217 <- bc_process(exp.name="20230217_mix0b")

t0.keep <- c("ALL-A","ALL-B","ALL-C",
             "DES1-B","DES1-D","DES1-E",
             "DES2-A","DES2-B","DES2-D",
             "EFF-A","EFF-B","EFF-C",
             "RPP-A","RPP-B","RPP-D",
             "minDES2-A","minDES2-B","minDES2-D")

mix.order <- c("ALL-A","ALL-B","ALL-C",
               "minDES2-A","minDES2-B","minDES2-D",
               "DES1-B","DES1-D","DES1-E",
               "DES2-A","DES2-B","DES2-D",
               "EFF-A","EFF-B","EFF-C",
               "RPP-A","RPP-B","RPP-D"
)

t0 <- rbind(EXP20230207, EXP20230217) %>%
  filter(mix_rep %in% t0.keep) %>%
  filter(drug_conc_timepoint != "none-none-20") %>%
  filter(drug_conc_timepoint != "TET-16-20")

t0$mix_rep <- factor(t0$mix_rep, levels=mix.order)

# Wrangle competitions
EXP20230806 <- bc_process(exp.name="20230806_allpools")
EXP20230806 <- EXP20230806 %>%
  subset(tech == "2") %>%
  subset(drug == "TIG" | drug == "none") %>%
  filter(mix_drug_conc_timepoint != "ALL-TIG-2-20")

EXP20230921A <- bc_process(exp.name="20230921A")
EXP20230921B <- bc_process(exp.name="20230921B")
EXP20230927 <- bc_process(exp.name="20230927")
EXP20230927.01 <- bc_process(exp.name="20230927-01")
EXP20230927.02 <- bc_process(exp.name="20230927-02") 

EXP20231128.01 <- bc_process(exp.name="20231128-01")
EXP20231128.02 <- bc_process(exp.name="20231128-02")

# Merge
merge <- rbind(t0, EXP20230921A, EXP20230921B, EXP20230806, EXP20230927.01, EXP20230927.02, EXP20231128.01, EXP20231128.02) 

# Remove those that needed to go for 48 h
all.merge <- merge %>%
  filter(mix_drug_conc_timepoint != "DES2-DOX-4-20") %>%
  filter(mix_drug_conc_timepoint != "DES2-MIN-4-20") %>%
  filter(mix_drug_conc_timepoint != "DES2-TET-4-20") %>%
  filter(mix_drug_conc_timepoint != "DES1-TIG-2-20") %>%
  filter(mix_drug_conc_timepoint != "DES1-TIG-0.5-20") %>%
  filter(mix_drug_conc_timepoint != "EFF-MIN-4-20") %>%
  filter(mix_drug_conc_timepoint != "RPP-MIN-4-20")


# # Add new column for plotting
# all.merge <- all.merge %>% mutate(plotlabel = plot.name[as.character(drug_conc_timepoint)])


# ### Barcode counts -------------------------------------------------------------
# p00.bc.count <- ggplot(EXP20231128.02, aes(x=mix_drug_conc_timepoint_rep, y=count)) + #all.merge
#   geom_bar(stat="identity") +
#   theme_plot() +
#   geom_hline(yintercept=10000)


### ALL competitions -----------------------------------------------------------

# Data wrangling 
ALL <- all.merge %>%
  filter(mix == "ALL") %>%
  filter(drug_conc_timepoint != "none-none-20")

# Plots
p01.ALL.TET.EFF <- comp_plot(mix.show=ALL, drug.show="TET", mech.show="EFF")
p01.ALL.DOX.RPP <- comp_plot(mix.show=ALL, drug.show="DOX", mech.show="RPP")
p01.ALL.MIN.RPP <- comp_plot(mix.show=ALL, drug.show="MIN", mech.show="RPP")
p01.ALL.TIG.DES1 <- comp_plot(mix.show=ALL, drug.show="TIG", mech.show="DES1")

# # Export plots 
# 
# # ALL + TET
# export_plot(p01.ALL.TET.EFF, fname="ALL-TET_v2", w=3, h=2, fmt="png")
# export_plot(p01.ALL.TET.EFF, fname="ALL-TET_v2", w=3, h=2, fmt="pdf")
# 
# # ALL + DOX
# export_plot(p01.ALL.DOX.RPP, fname="ALL-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p01.ALL.DOX.RPP, fname="ALL-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # ALL + MIN
# export_plot(p01.ALL.MIN.RPP, fname="ALL-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p01.ALL.MIN.RPP, fname="ALL-MIN_v2", w=3, h=2, fmt="pdf")
# 
# # ALL + TIG
# export_plot(p01.ALL.TIG.DES1, fname="ALL-TIG_v2", w=3, h=2, fmt="png")
# export_plot(p01.ALL.TIG.DES1, fname="ALL-TIG_v2", w=3, h=2, fmt="pdf")

# ALL quantification 
ALL.TET.EFF.df <- comp_analyze(mix.show=ALL, drug.show="TET", mech.show="EFF")
ALL.DOX.RPP.df <- comp_analyze(mix.show=ALL, drug.show="DOX", mech.show="RPP")
ALL.MIN.RPP.df <- comp_analyze(mix.show=ALL, drug.show="MIN", mech.show="RPP")
ALL.TIG.DES1.df <- comp_analyze(mix.show=ALL, drug.show="TIG", mech.show="DES1")


### Relab MIC comparison -------------------------------------------------------

# TET
ALL.TET.df <- comp_analyze(mix.show=ALL, drug.show="TET")
ALL.TET.merge <- merge(ALL.TET.df[[5]], MIC.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p10.TET.MIC.relab <- ggplot(ALL.TET.merge, aes(x=TET, y=mean_gene_frac, color=gene, group=drug_conc_timepoint)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.TET.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2, size=0.5) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel abundance", x="TET MIC foldchange")

print(p10.TET.MIC.relab)

# DOX
ALL.DOX.df <- comp_analyze(mix.show=ALL, drug.show="DOX")
ALL.DOX.merge <- merge(ALL.DOX.df[[5]], MIC.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p10.DOX.MIC.relab <- ggplot(ALL.DOX.merge, aes(x=DOX, y=mean_gene_frac, color=gene, group=drug_conc_timepoint)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.DOX.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2, size=0.5) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel. abundance", x="DOX MIC foldchange")

# MIN
ALL.MIN.df <- comp_analyze(mix.show=ALL, drug.show="MIN")
ALL.MIN.merge <- merge(ALL.MIN.df[[5]], MIC.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p10.MIN.MIC.relab <- ggplot(ALL.MIN.merge, aes(x=MIN, y=mean_gene_frac, color=gene, group=drug_conc_timepoint)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.MIN.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2, size=0.5) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel. abundance", x="MIN MIC foldchange")

# TIG
ALL.TIG.df <- comp_analyze(mix.show=ALL, drug.show="TIG")
ALL.TIG.merge <- merge(ALL.TIG.df[[5]], MIC.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p10.TIG.MIC.relab <- ggplot(ALL.TIG.merge, aes(x=TIG, y=mean_gene_frac, color=gene, group=drug_conc_timepoint)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.TIG.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2, size=0.5) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel. abundance", x="TIG MIC foldchange")

p10.MIC.relab.merge <- ggarrange(p10.TET.MIC.relab, p10.DOX.MIC.relab, p10.MIN.MIC.relab, p10.TIG.MIC.relab)

# export_plot(p10.MIC.relab.merge, fname="ALL-relab-MIC_v4", w=8, h=3, fmt="png")
# export_plot(p10.MIC.relab.merge, fname="ALL-relab-MIC_v4", w=8, h=3, fmt="pdf")


### Relab GR comparison --------------------------------------------------------

# TET
ALL.TET.df <- comp_analyze(mix.show=ALL, drug.show="TET")
ALL.TET.GR.merge <- merge(ALL.TET.df[[5]], GR.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p20.TET.GR.relab <- ggplot(ALL.TET.GR.merge, aes(x=Average, y=mean_gene_frac, color=gene)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.TET.GR.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.01, size=0.5) + #
  geom_errorbarh(ALL.TET.GR.merge, mapping=aes(xmax=Average+StDev, xmin=Average-StDev), height=0.01) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel abundance", x="TET MIC foldchange")

# DOX
ALL.DOX.df <- comp_analyze(mix.show=ALL, drug.show="DOX")
ALL.DOX.GR.merge <- merge(ALL.DOX.df[[5]], GR.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p20.DOX.GR.relab <- ggplot(ALL.DOX.GR.merge, aes(x=Average, y=mean_gene_frac, color=gene)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.DOX.GR.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.01, size=0.5) + #
  geom_errorbarh(ALL.DOX.GR.merge, mapping=aes(xmax=Average+StDev, xmin=Average-StDev), height=0.01) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel abundance", x="DOX MIC foldchange")

# MIN
ALL.MIN.df <- comp_analyze(mix.show=ALL, drug.show="MIN")
ALL.MIN.GR.merge <- merge(ALL.MIN.df[[5]], GR.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p20.MIN.GR.relab <- ggplot(ALL.MIN.GR.merge, aes(x=Average, y=mean_gene_frac, color=gene)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.MIN.GR.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.01, size=0.5) + #
  geom_errorbarh(ALL.MIN.GR.merge, mapping=aes(xmax=Average+StDev, xmin=Average-StDev), height=0.01) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel abundance", x="MIN MIC foldchange")

# TIG
ALL.TIG.df <- comp_analyze(mix.show=ALL, drug.show="TIG")
ALL.TIG.GR.merge <- merge(ALL.TIG.df[[5]], GR.results, by=c("gene")) %>%
  filter(drug_conc_timepoint != "none-none-0")

p20.TIG.GR.relab <- ggplot(ALL.TIG.GR.merge, aes(x=Average, y=mean_gene_frac, color=gene)) +
  #geom_smooth(aes(group=drug_conc_timepoint), method="lm") +
  geom_point() +
  theme_plot() +
  theme(
    axis.text.x = element_text(size=6, angle = 0, vjust=1, hjust=0.5), #, hjust=1, vjust=0.5
  ) +
  geom_errorbar(ALL.TIG.GR.merge, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.01, size=0.5) + #
  geom_errorbarh(ALL.TIG.GR.merge, mapping=aes(xmax=Average+StDev, xmin=Average-StDev), height=0.01) +
  scale_color_manual(values=gene.pal) +
  facet_grid(~drug_conc_timepoint) +
  labs(y="rel abundance", x="TIG MIC foldchange")

p20.GR.relab.merge <- ggarrange(p20.TET.GR.relab, p20.DOX.GR.relab, p20.MIN.GR.relab, p20.TIG.GR.relab)
print(p20.GR.relab.merge)

# export_plot(p20.TET.GR.relab, fname="ALL-relab-GR_v4", w=6, h=2, fmt="png")
# export_plot(p20.TET.GR.relab, fname="ALL-relab-GR_v4", w=6, h=2, fmt="pdf")


### minDES2 competitions -------------------------------------------------------

# minDES2 data wrangling 
minDES2 <- all.merge %>%
  filter(mix == "minDES2") %>%
  filter(drug_conc_timepoint != "none-none-20") %>%
  filter(drug_conc_timepoint != "TIG-0.5-48") %>% #0.5-20
  filter(drug_conc_timepoint != "TIG-2-20") 


# minDES2 plots 
p01.minDES2.TET.EFF <- comp_plot(mix.show=minDES2, drug.show="TET", mech.show="EFF")
p01.minDES2.DOX.RPP <- comp_plot(mix.show=minDES2, drug.show="DOX", mech.show="RPP")
p01.minDES2.MIN.RPP <- comp_plot(mix.show=minDES2, drug.show="MIN", mech.show="RPP")
p01.minDES2.TIG.DES1 <- comp_plot(mix.show=minDES2, drug.show="TIG", mech.show="DES1", xaxis.order=minDES2.TIG.order)

# # Export plots 
# 
# # minDES2 + TET
# export_plot(p01.minDES2.TET.EFF, fname="minDES2-TET_v2", w=3, h=2, fmt="png")
# export_plot(p01.minDES2.TET.EFF, fname="minDES2-TET_v2", w=3, h=2, fmt="pdf")
# 
# # minDES2 + DOX
# export_plot(p01.minDES2.DOX.RPP, fname="minDES2-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p01.minDES2.DOX.RPP, fname="minDES2-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # minDES2 + MIN
# export_plot(p01.minDES2.MIN.RPP, fname="minDES2-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p01.minDES2.MIN.RPP, fname="minDES2-MIN_v2", w=3, h=2, fmt="pdf")
# 
# # minDES2 + TIG
# export_plot(p01.minDES2.TIG.DES1, fname="minDES2-TIG_v2", w=3, h=2, fmt="png")
# export_plot(p01.minDES2.TIG.DES1, fname="minDES2-TIG_v2", w=3, h=2, fmt="pdf")

# minDES2 quantification 
minDES2.TET.EFF.df <- comp_analyze(mix.show=minDES2, drug.show="TET", mech.show="EFF")
minDES2.DOX.RPP.df <- comp_analyze(mix.show=minDES2, drug.show="DOX", mech.show="RPP")
minDES2.MIN.RPP.df <- comp_analyze(mix.show=minDES2, drug.show="MIN", mech.show="RPP")
minDES2.TIG.DES1.df <- comp_analyze(mix.show=minDES2, drug.show="TIG", mech.show="DES1")


### Mechanism-specific competitions --------------------------------------------

# DES1 data wrangling
DES1 <- all.merge %>%
  filter(mix == "DES1") %>%
  filter(drug_conc_timepoint != "none-none-20")

DES1$count[DES1$count <= 50] <- 0

# DES1 plots
p02.DES1.TET <- comp_plot_mech(mix.show=DES1, drug.show="TET", mech.show="DES1")
p02.DES1.DOX <- comp_plot_mech(mix.show=DES1, drug.show="DOX", mech.show="DES1")
p02.DES1.MIN <- comp_plot_mech(mix.show=DES1, drug.show="MIN", mech.show="DES1")
p02.DES1.TIG <- comp_plot_mech(mix.show=DES1, drug.show="TIG", mech.show="DES1", xaxis.order=DES1.TIG.order)

# # Export plots
# # DES1 + TET
# export_plot(p02.DES1.TET, fname="DES1-TET_v2", w=3, h=2, fmt="png")
# export_plot(p02.DES1.TET, fname="DES1-TET_v2", w=3, h=2, fmt="pdf")
# 
# # DES1 + DOX
# export_plot(p02.DES1.DOX, fname="DES1-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p02.DES1.DOX, fname="DES1-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # DES1 + MIN
# export_plot(p02.DES1.MIN, fname="DES1-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p02.DES1.MIN, fname="DES1-MIN_v2", w=3, h=2, fmt="pdf")
# 
# # DES1 + TIG
# export_plot(p02.DES1.TIG, fname="DES1-TIG_v2", w=3, h=2, fmt="png")
# export_plot(p02.DES1.TIG, fname="DES1-TIG_v2", w=3, h=2, fmt="pdf")


# DES1 quantification
DES1.TET.df <- comp_analyze(mix.show=DES1, drug.show="TET", mech.show="DES1")
DES1.DOX.df <- comp_analyze(mix.show=DES1, drug.show="DOX", mech.show="DES1")
DES1.MIN.df <- comp_analyze(mix.show=DES1, drug.show="MIN", mech.show="DES1")
DES1.TIG.df <- comp_analyze(mix.show=DES1, drug.show="TIG", mech.show="DES1")

### EFF competitions -----------------------------------------------------------

# EFF data wrangling
EFF <- all.merge %>%
  filter(mix == "EFF") %>%
  filter(drug_conc_timepoint != "none-none-20")

EFF$count[EFF$count <= 50] <- 0

# EFF plots
p03.EFF.TET <- comp_plot_mech(mix.show=EFF, drug.show="TET", mech.show="EFF")
p03.EFF.DOX <- comp_plot_mech(mix.show=EFF, drug.show="DOX", mech.show="EFF")
p03.EFF.MIN <- comp_plot_mech(mix.show=EFF, drug.show="MIN", mech.show="EFF", xaxis.order=EFF.MIN.order)

# # Export plots
# # EFF + TET
# export_plot(p03.EFF.TET, fname="EFF-TET_v2", w=3, h=2, fmt="png")
# export_plot(p03.EFF.TET, fname="EFF-TET_v2", w=3, h=2, fmt="pdf")
# 
# # EFF + DOX
# export_plot(p03.EFF.DOX, fname="EFF-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p03.EFF.DOX, fname="EFF-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # EFF + MIN
# export_plot(p03.EFF.MIN, fname="EFF-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p03.EFF.MIN, fname="EFF-MIN_v2", w=3, h=2, fmt="pdf")


# EFF quantification
EFF.TET.df <- comp_analyze(mix.show=EFF, drug.show="TET", mech.show="EFF")
EFF.DOX.df <- comp_analyze(mix.show=EFF, drug.show="DOX", mech.show="EFF")
EFF.MIN.df <- comp_analyze(mix.show=EFF, drug.show="MIN", mech.show="EFF")


### DES2 competitions ----------------------------------------------------------

# DES2 data wrangling
DES2 <- all.merge %>%
  filter(mix == "DES2") %>%
  filter(drug_conc_timepoint != "none-none-20")

DES2$count[DES2$count <= 60] <- 0

# DES2 plots
p04.DES2.TET <- comp_plot_mech(mix.show=DES2, drug.show="TET", mech.show="DES2", xaxis.order=DES2.TET.order)
p04.DES2.DOX <- comp_plot_mech(mix.show=DES2, drug.show="DOX", mech.show="DES2", xaxis.order=DES2.DOX.order)
p04.DES2.MIN <- comp_plot_mech(mix.show=DES2, drug.show="MIN", mech.show="DES2", xaxis.order=DES2.MIN.order)

# # Export plots
# # DES2 + TET
# export_plot(p04.DES2.TET, fname="DES2-TET_v2", w=3, h=2, fmt="png")
# export_plot(p04.DES2.TET, fname="DES2-TET_v2", w=3, h=2, fmt="pdf")
# 
# # DES2 + DOX
# export_plot(p04.DES2.DOX, fname="DES2-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p04.DES2.DOX, fname="DES2-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # DES2 + MIN
# export_plot(p04.DES2.MIN, fname="DES2-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p04.DES2.MIN, fname="DES2-MIN_v2", w=3, h=2, fmt="pdf")


# DES2 quantification
DES2.TET.df <- comp_analyze(mix.show=DES2, drug.show="TET", mech.show="DES2")
DES2.DOX.df <- comp_analyze(mix.show=DES2, drug.show="DOX", mech.show="DES2")
DES2.MIN.df <- comp_analyze(mix.show=DES2, drug.show="MIN", mech.show="DES2")


### RPP competitions -----------------------------------------------------------

# RPP data wrangling
RPP <- all.merge %>%
  filter(mix == "RPP") %>%
  filter(drug_conc_timepoint != "none-none-20")

RPP$count[RPP$count <= 50] <- 0

# RPP plots
p05.RPP.TET <- comp_plot_mech(mix.show=RPP, drug.show="TET", mech.show="RPP")
p05.RPP.DOX <- comp_plot_mech(mix.show=RPP, drug.show="DOX", mech.show="RPP")
p05.RPP.MIN <- comp_plot_mech(mix.show=RPP, drug.show="MIN", mech.show="RPP", xaxis.order=RPP.MIN.order)

# # Export plots
# # RPP + TET
# export_plot(p05.RPP.TET, fname="RPP-TET_v2", w=3, h=2, fmt="png")
# export_plot(p05.RPP.TET, fname="RPP-TET_v2", w=3, h=2, fmt="pdf")
# 
# # RPP + DOX
# export_plot(p05.RPP.DOX, fname="RPP-DOX_v2", w=3, h=2, fmt="png")
# export_plot(p05.RPP.DOX, fname="RPP-DOX_v2", w=3, h=2, fmt="pdf")
# 
# # RPP + MIN
# export_plot(p05.RPP.MIN, fname="RPP-MIN_v2", w=3, h=2, fmt="png")
# export_plot(p05.RPP.MIN, fname="RPP-MIN_v2", w=3, h=2, fmt="pdf")


# RPP quantification
RPP.TET.df <- comp_analyze(mix.show=RPP, drug.show="TET", mech.show="RPP")
RPP.DOX.df <- comp_analyze(mix.show=RPP, drug.show="DOX", mech.show="RPP")
RPP.MIN.df <- comp_analyze(mix.show=RPP, drug.show="MIN", mech.show="RPP")


### Antibiotic-free competitions (none) ----------------------------------------

ALL.none <- all.merge %>%
  filter(mix == "ALL") %>%
  filter(drug == "none")

p00a.ALL.none.barplot <- ggplot(ALL.none, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.none.change <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.none.relab <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p00b.ALL.none.mech.per <- ggplot(ALL.none.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.none.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

ALL.none.relab.gene <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
  summarize(tot_gene = sum(count)) %>%
  mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
  ungroup() %>%
  group_by(drug_conc_timepoint, gene) %>% # Get summary stats of all reps combined
  summarize(mean_gene_frac = mean(gene_frac), sd=sd(gene_frac))

p00.ALL.none.merge <- ggarrange(p00a.ALL.none.barplot, p00b.ALL.none.mech.per, ncol=3)

# fname.ALL.none.pdf <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none_v1.pdf")
# fname.ALL.none.png <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none_v1.png")
# ggsave(p00.ALL.none.merge, filename = fname.ALL.none.pdf , device=cairo_pdf, height=2.05, width=1.35) #cairo_pdf
#ggsave(p00.ALL.none.merge, filename = fname.ALL.none.png , device=png, height=1.75, width=2) #cairo_pdf

### More plots for none comparisons --------------------------------------------

# Specific comparisons
ALL.none.relab.gene.stats <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
  summarize(tot_gene = sum(count)) %>%
  mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
  ungroup() %>%
  group_by(gene) %>%
  pairwise_t_test(gene_frac ~ drug_conc_timepoint, p.adjust.method = "bonferroni")
ALL.none.relab.gene.stats.sub <- ALL.none.relab.gene.stats %>%
  subset(p.adj.signif != "ns")

# Add to get indiv points for barplots
ALL.none.relab.gene.rep <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, gene) %>%  # Get frac ea gene per ea rep
  summarize(tot_gene = sum(count)) %>%
  mutate(gene_frac = tot_gene/sum(tot_gene)) %>% # Subtract ea gene_frac by the gene_frac of t0 (== gene_frac[1])
  ungroup()

p00c.ALL.none.gene.bar <- ggplot(ALL.none.relab.gene, aes(x=drug_conc_timepoint, y=mean_gene_frac, fill=gene)) +
  geom_bar(stat="identity") +
  geom_errorbar(ALL.none.relab.gene, mapping=aes(ymin=mean_gene_frac-sd, ymax=mean_gene_frac+sd), width=0.2) +
  geom_jitter(data=ALL.none.relab.gene.rep, aes(x=drug_conc_timepoint, y=gene_frac), alpha=0.4, width=0.3) +
  theme_plot() +
  stat_pvalue_manual(ALL.none.relab.gene.stats.sub, label="p.adj.signif", y.position=0.08) +
  ylim(0, 0.1) +
  facet_wrap(~ gene) +
  scale_fill_manual(values=gene.pal) +
  labs(y="avg. abundance")

print(p00c.ALL.none.gene.bar)

# fname.ALL.none.gene.pdf <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none-gene_S03B_v1.pdf")
# fname.ALL.none.gene.png <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none-gene_S03B_v1.png")
# ggsave(p00c.ALL.none.gene.bar, filename = fname.ALL.none.gene.pdf , device=cairo_pdf) #cairo_pdf
# ggsave(p00c.ALL.none.gene.bar, filename = fname.ALL.none.gene.png , device=png) #cairo_pdf

ALL.none.relab.stats <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(mech) %>%
  pairwise_t_test(mech_frac ~ drug_conc_timepoint, p.adjust.method = "bonferroni")
ALL.none.relab.mech.stats.sub <- ALL.none.relab.stats %>%
  subset(p.adj.signif != "ns")

# Add to get indiv points for barplot
ALL.none.relab.rep <- ALL.none %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup()

p00d.ALL.none.mech.bar <- ggplot(ALL.none.relab, aes(x=drug_conc_timepoint, y=mean_mech_frac, fill=mech)) +
  geom_bar(stat="identity") +
  geom_errorbar(ALL.none.relab, mapping=aes(ymin=mean_mech_frac-sd, ymax=mean_mech_frac+sd), width=0.2) +
  geom_jitter(data=ALL.none.relab.rep, aes(x=drug_conc_timepoint, y=mech_frac), alpha=0.4, width=0.3) +
  theme_plot() +
  stat_pvalue_manual(ALL.none.relab.mech.stats.sub, label="p.adj.signif", y.position=0.4) +
  ylim(0, 0.5) +
  facet_wrap(~ mech) +
  scale_fill_manual(values=mech.pal) +
  labs(y="avg. abundance")

print(p00d.ALL.none.mech.bar)

# fname.ALL.none.mech.pdf <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none-mech_S03A_v1.pdf")
# fname.ALL.none.mech.png <- stringr::str_interp("${out.path}/${fln.prfx}-ALL-none-mech_S03A_v1.png")
# ggsave(p00d.ALL.none.mech.bar, filename = fname.ALL.none.mech.pdf , device=cairo_pdf, width=, height=2) #cairo_pdf
# ggsave(p00d.ALL.none.mech.bar, filename = fname.ALL.none.mech.png , device=png, width=2, height=2) #cairo_pdf


# ### ATC experiments ------------------------------------------------------------
# ALL <- all.merge %>%
#   filter(mix == "ALL")
# 
# p99.ALL.ATC.DES1 <- comp_plot(mix.show=ALL, drug.show="ATC", mech.show="DES1")
# 
# print(p99.ALL.ATC.DES1)
# 
# # ALL + TET
# export_plot(p99.ALL.ATC.DES1, fname="ALL-ATC", w=3, h=2, fmt="png")
# export_plot(p99.ALL.ATC.DES1, fname="ALL-ATC", w=3, h=2, fmt="pdf")
# 
# # ALL quantification 
# ALL.ATC.DES1.df <- comp_analyze(mix.show=ALL, drug.show="ATC", mech.show="DES1")
# 
# 
# # Breakdown by tech rep 
# EXP20231128_01_atc <- bc_process(exp.name="20231128-01") %>%
#   subset(drug == "ATC")
# 
# p99.ALL.ATC.reps <- ggplot(EXP20231128_01_atc, aes(x=mix_rep, y=count, fill=gene)) +
#   geom_bar(stat="identity", position="fill") +
#   theme_plot() +
#   scale_fill_manual(values=gene.pal) +
#   facet_grid(cols=vars(drug_conc)) +
#   labs(y="relative abundance")
# 
# # export_plot(p99.ALL.ATC.reps, fname="ALL-ATC-replicates", w=3, h=2, fmt="png")
# # export_plot(p99.ALL.ATC.reps, fname="ALL-ATC-replicates", w=3, h=2, fmt="pdf")


### Transfer experiments -------------------------------------------------------

ALL.transf <- all.merge %>%
  filter(mix == "ALL") %>%
  filter(drug == "transf" | drug =="none") %>%
  filter(drug_conc_timepoint != "none-none-20")

### TET start 
ALL.TET.transf.01 <- ALL.transf %>%
  filter(drug_conc == "transf-TET" | drug_conc == "transf-TET-DOX" | drug_conc == "transf-TET-DOX-TIG" | drug_conc == "none-none")

p50a.transf.TET.01.barplot <- ggplot(ALL.TET.transf.01, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.TET.transf.01.change <- ALL.TET.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.TET.transf.01.relab <- ALL.TET.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p50b.transf.TET.01.mech.per <- ggplot(ALL.TET.transf.01.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.TET.transf.01.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p50.TET.transf.01.merge <- ggarrange(p50a.transf.TET.01.barplot, p50b.transf.TET.01.mech.per, ncol=3)

ALL.TET.transf.02 <- ALL.transf %>%
  filter(drug_conc == "transf-TET" | drug_conc == "transf-TET-TIG" | drug_conc == "transf-TET-TIG-DOX" | drug_conc == "none-none")

p51a.transf.TET.02.barplot <- ggplot(ALL.TET.transf.02, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.TET.transf.02.change <- ALL.TET.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.TET.transf.02.relab <- ALL.TET.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p51b.transf.TET.02.mech.per <- ggplot(ALL.TET.transf.02.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.TET.transf.02.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p51.TET.transf.02.merge <- ggarrange(p51a.transf.TET.02.barplot, p51b.transf.TET.02.mech.per, ncol=3)

### DOX start 
ALL.DOX.transf.01 <- ALL.transf %>%
  filter(drug_conc == "transf-DOX" | drug_conc == "transf-DOX-TET" | drug_conc == "transf-DOX-TET-TIG" | drug_conc == "none-none")

p52a.transf.DOX.01.barplot <- ggplot(ALL.DOX.transf.01, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.DOX.transf.01.change <- ALL.DOX.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.DOX.transf.01.relab <- ALL.DOX.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p52b.transf.DOX.01.mech.per <- ggplot(ALL.DOX.transf.01.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.DOX.transf.01.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p52.DOX.transf.01.merge <- ggarrange(p52a.transf.DOX.01.barplot, p52b.transf.DOX.01.mech.per , ncol=3)

ALL.DOX.transf.02 <- ALL.transf %>%
  filter(drug_conc == "transf-DOX" | drug_conc == "transf-DOX-TIG" | drug_conc == "transf-DOX-TIG-TET" | drug_conc == "none-none")

p53a.transf.DOX.02.barplot <- ggplot(ALL.DOX.transf.02, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.DOX.transf.02.change <- ALL.DOX.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.DOX.transf.02.relab <- ALL.DOX.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p53b.transf.DOX.02.mech.per <- ggplot(ALL.DOX.transf.02.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.DOX.transf.02.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p53.DOX.transf.02.merge <- ggarrange(p53a.transf.DOX.02.barplot, p53b.transf.DOX.02.mech.per, ncol=3)


### TIG start 
ALL.TIG.transf.01 <- ALL.transf %>%
  filter(drug_conc == "transf-TIG" | drug_conc == "transf-TIG-TET" | drug_conc == "transf-TIG-TET-DOX" | drug_conc == "none-none" )

p54a.transf.TIG.01.barplot <- ggplot(ALL.TIG.transf.01, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.TIG.transf.01.change <- ALL.TIG.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.TIG.transf.01.relab <- ALL.TIG.transf.01 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p54b.transf.TIG.01.mech.per <- ggplot(ALL.TIG.transf.01.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.TIG.transf.01.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p54.TIG.transf.01.merge <- ggarrange(p54a.transf.TIG.01.barplot, p54b.transf.TIG.01.mech.per , ncol=3)

ALL.TIG.transf.02 <- ALL.transf %>%
  filter(drug_conc == "transf-TIG" | drug_conc == "transf-TIG-DOX" | drug_conc == "transf-TIG-DOX-TET" | drug_conc == "none-none")

p55a.transf.TIG.02.barplot <- ggplot(ALL.TIG.transf.02, aes(x=drug_conc_timepoint, y=count, fill=gene)) +
  geom_bar(stat="identity", position="fill") +
  theme_plot() +
  #theme(legend.position="right") +
  scale_fill_manual(values=gene.pal) +
  labs(y="rel. abundance")

# Per change in mech rel ab
ALL.TIG.transf.02.change <- ALL.TIG.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  group_by(rep, mech) %>% # Subtract ea mech_frac by the mech_frac of t0 (== mech_frac[1])
  mutate(change_per = ((mech_frac - mech_frac[1]) / mech_frac[1]) * 100) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>% # Get summary stats of all reps combined
  summarize(mean_mech_per = mean(change_per), sd=sd(change_per)) %>%
  ungroup()

ALL.TIG.transf.02.relab <- ALL.TIG.transf.02 %>%
  group_by(drug_conc_timepoint, rep, mech) %>% # Get frac ea mech per ea rep
  summarize(tot_mech = sum(count)) %>%
  mutate(mech_frac = tot_mech/sum(tot_mech)) %>%
  ungroup() %>%
  group_by(drug_conc_timepoint, mech) %>%
  summarize(mean_mech_frac = mean(mech_frac), sd=sd(mech_frac))

p55b.transf.TIG.02.mech.per <- ggplot(ALL.TIG.transf.02.change, aes(x=drug_conc_timepoint, y=mean_mech_per, color=mech, group=mech)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(ALL.TIG.transf.02.change, mapping=aes(ymin=mean_mech_per-sd, ymax=mean_mech_per+sd), width=0.2, size=0.5) +
  theme_plot() +
  #theme(legend.position="right") +
  #ylim(-100, 350) +
  scale_color_manual(values=mech.pal) +
  labs(y="% change")

p55.TIG.transf.01.merge <- ggarrange(p55a.transf.TIG.02.barplot, p55b.transf.TIG.02.mech.per , ncol=3)

