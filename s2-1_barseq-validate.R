
###---Load data-----------------------------------------------------------------

#B22 == pipetting error duplicated in tech rep rxns 003 and 004
val02_10k <- rare_process("20230213_val02_rarefaction/20230213_val02-10k")
  val02_10k_sub <- filter(val02_10k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_10k_sub <- lm(relab ~ measab, data = val02_10k_sub)
  val02_10k_sub2 <- val02_10k[!(val02_10k$barcode == "b22" & val02_10k$rxn == "PCR1-003") & !(val02_10k$barcode == "b22" & val02_10k$rxn == "PCR1-004") & !(val02_10k$rxn == "PCR1-005"),]
  val02_10k_sub3 <- val02_10k[!(val02_10k$barcode == "b22" & val02_10k$rxn == "PCR1-003") & !(val02_10k$rxn == "PCR1-004") & !(val02_10k$rxn == "PCR1-005"),]
val02_30k <- rare_process("20230213_val02_rarefaction/20230213_val02-30k")
  val02_30k_sub <- filter(val02_30k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_30k_sub <- lm(relab ~ measab, data = val02_30k_sub)
  val02_30k_sub2 <- val02_30k[!(val02_30k$barcode == "b22" & val02_30k$rxn == "PCR1-003") & !(val02_30k$barcode == "b22" & val02_30k$rxn == "PCR1-004") & !(val02_30k$rxn == "PCR1-005"),]
  val02_30k_sub3 <- val02_30k[!(val02_30k$barcode == "b22" & val02_30k$rxn == "PCR1-003") & !(val02_30k$rxn == "PCR1-004") & !(val02_30k$rxn == "PCR1-005"),]
val02_50k <- rare_process("20230213_val02_rarefaction/20230213_val02-50k")
  val02_50k_sub <- filter(val02_50k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_50k_sub <- lm(relab ~ measab, data = val02_50k_sub)
  val02_50k_sub2 <- val02_50k[!(val02_50k$barcode == "b22" & val02_50k$rxn == "PCR1-003") & !(val02_50k$barcode == "b22" & val02_50k$rxn == "PCR1-004") & !(val02_50k$rxn == "PCR1-005"),]
  val02_50k_sub3 <- val02_50k[!(val02_50k$barcode == "b22" & val02_50k$rxn == "PCR1-003") & !(val02_50k$rxn == "PCR1-004") & !(val02_50k$rxn == "PCR1-005"),]
val02_100k <- rare_process("20230213_val02_rarefaction/20230213_val02-100k")
  val02_100k_sub <- filter(val02_100k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_100k_sub <- lm(relab ~ measab, data = val02_100k_sub)
  val02_100k_sub2 <- val02_100k[!(val02_100k$barcode == "b22" & val02_100k$rxn == "PCR1-003") & !(val02_100k$barcode == "b22" & val02_100k$rxn == "PCR1-004") & !(val02_100k$rxn == "PCR1-005"),]
  val02_100k_sub3 <- val02_100k[!(val02_100k$barcode == "b22" & val02_100k$rxn == "PCR1-003") & !(val02_100k$rxn == "PCR1-004") & !(val02_100k$rxn == "PCR1-005"),]
val02_150k <- rare_process("20230213_val02_rarefaction/20230213_val02-150k")
  val02_150k_sub <- filter(val02_150k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_150k_sub <- lm(relab ~ measab, data = val02_150k_sub)
  val02_150k_sub2 <- val02_150k[!(val02_150k$barcode == "b22" & val02_150k$rxn == "PCR1-003") & !(val02_150k$barcode == "b22" & val02_150k$rxn == "PCR1-004") & !(val02_150k$rxn == "PCR1-005"),]
  val02_150k_sub3 <- val02_150k[!(val02_150k$barcode == "b22" & val02_150k$rxn == "PCR1-003") & !(val02_150k$rxn == "PCR1-004") & !(val02_150k$rxn == "PCR1-005"),]
val02_200k <- rare_process("20230213_val02_rarefaction/20230213_val02-200k")
  val02_200k_sub <- filter(val02_200k, rxn != "PCR1-003" & rxn != "PCR1-004" & rxn != "PCR1-005")
  fit_200k_sub <- lm(relab ~ measab, data = val02_200k_sub)
  val02_200k_sub2 <- val02_200k[!(val02_200k$barcode == "b22" & val02_200k$rxn == "PCR1-003") & !(val02_200k$barcode == "b22" & val02_200k$rxn == "PCR1-004") & !(val02_200k$rxn == "PCR1-005"),]
  val02_200k_sub3 <- val02_200k[!(val02_200k$barcode == "b22" & val02_200k$rxn == "PCR1-003") & !(val02_200k$rxn == "PCR1-004") & !(val02_200k$rxn == "PCR1-005"),]

# ###---JUST 001 & 002----------------------------------------------------------
# p01.val02.10k <- ggplot(val02_10k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "10k reads")
#   
# print(p01.val02.10k)
# 
# p01.val02.30k <- ggplot(val02_30k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "30k reads")
#   
# p01.val02.50k <- ggplot(val02_50k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "50k reads")
#   
# p01.val02.100k <- ggplot(val02_100k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "100k reads")
#   
#   
# p01.val02.150k <- ggplot(val02_150k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "150k reads")
# 
# p01.val02.200k <- ggplot(val02_200k_sub, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "200k reads")
# 
# val02.rare.fig <- ggarrange(p01.val02.10k, p01.val02.30k, p01.val02.50k, p01.val02.100k, p01.val02.150k, p01.val02.200k, ncol=3, nrow=2)
# 
# print(val02.rare.fig)
# 
# fname.val02.rare.fig <- stringr::str_interp("${out.path}/${fln.prfx}-val02_rare_v1.pdf")
# ggsave(val02.rare.fig, filename = fname.val02.rare.fig, device=cairo_pdf, width=7) #width=5 

###---MASK B22s-----------------------------------------------------------------

p02.val02.10k <- ggplot(val02_10k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) + #label.y = 0.01, 
  labs(title = "10k reads")

print(p02.val02.10k)

p02.val02.30k <- ggplot(val02_30k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "30k reads")

p02.val02.50k <- ggplot(val02_50k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),        
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) + #0.1
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "50k reads")

p02.val02.100k <- ggplot(val02_100k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),        
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "100k reads")


p02.val02.150k <- ggplot(val02_150k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),        
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "150k reads")

p02.val02.200k <- ggplot(val02_200k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.75) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),        
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  scale_y_continuous(trans='log10', limits=c(0.0001, 0.3)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "200k reads")

val02.rare.p02 <- ggarrange(p02.val02.10k, p02.val02.30k, p02.val02.50k, p02.val02.100k, p02.val02.150k, p02.val02.200k, ncol=3, nrow=2)

p02.val02.50k.plot <- ggplot(val02_50k_sub2, aes(x=relab, y=measab)) +
  geom_abline(intercept=0, slope=1, color="gray70") +
  geom_point(size=0.8) +
  stat_smooth(method = "lm") +
  theme_pub() +
  theme(
    aspect.ratio=1,
    axis.line.x =element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),        
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "right",        
    axis.ticks = element_line(color = "black", size=0.75)
  ) +
  scale_x_continuous(trans='log10', limits=c(0.0001, 1)) + #0.1
  scale_y_continuous(trans='log10', limits=c(0.0001, 1)) +
  stat_regline_equation(aes(label = ..eq.label..), size=2) +
  stat_regline_equation(aes(label = ..rr.label..), size=2) +
  labs(title = "50k reads")

# fname.val02.p02 <- stringr::str_interp("${out.path}/${fln.prfx}-val02_rare_v3.pdf")
# ggsave(val02.rare.p02, filename = fname.val02.p02, device=cairo_pdf, width=4, height=3) #width=5 

# ###---KEEP ALL------------------------------------------------------------------
# p03.val02.10k <- ggplot(val02_10k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "10k reads")
#   
# print(p03.val02.10k)
# 
# p03.val02.30k <- ggplot(val02_30k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "30k reads")
#   
# p03.val02.50k <- ggplot(val02_50k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "50k reads")
#   
# p03.val02.100k <- ggplot(val02_100k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "100k reads")
#   
#   
# p03.val02.150k <- ggplot(val02_150k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "150k reads")
# 
# p03.val02.200k <- ggplot(val02_200k, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "200k reads")
# 
# val02.rare.p03 <- ggarrange(p03.val02.10k, p03.val02.30k, p03.val02.50k, p03.val02.100k, p03.val02.150k, p03.val02.200k, ncol=3, nrow=2)
# 
# print(val02.rare.p03)
# 
# fname.val02.p03 <- stringr::str_interp("${out.path}/${fln.prfx}-val02_rare_v3.pdf")
# ggsave(val02.rare.p03, filename = fname.val02.p03, device=cairo_pdf, width=7) #width=5 
# 
# 
# ###---Just 001-003 AND MASKING B22----------------------------------------------
# p04.val02.10k <- ggplot(val02_10k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "10k reads")
#   
# print(p04.val02.10k)
# 
# p04.val02.30k <- ggplot(val02_30k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "30k reads")
#   
# p04.val02.50k <- ggplot(val02_50k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "50k reads")
#   
# p04.val02.100k <- ggplot(val02_100k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "100k reads")
#   
#   
# p04.val02.150k <- ggplot(val02_150k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "150k reads")
# 
# p04.val02.200k <- ggplot(val02_200k_sub3, aes(x=relab, y=measab)) +
#   geom_abline(intercept=0, slope=1, color="gray70") +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   # xlim(0, 0.3) +
#   # ylim(0, 0.3) +
#   theme_pub() +
#   theme(
#         panel.border = element_rect(colour = "black", fill=NA, size=1),
#         legend.position = "right",        
#         axis.ticks = element_line(color = "black", size=0.75)
#         ) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.002, aes(label = ..rr.label..)) +
#   labs(title = "200k reads")
# 
# val02.rare.p04 <- ggarrange(p04.val02.10k, p04.val02.30k, p04.val02.50k, p04.val02.100k, p04.val02.150k, p04.val02.200k, ncol=3, nrow=2)
# 
# print(val02.rare.p04)
# 
# fname.val02.p04 <- stringr::str_interp("${out.path}/${fln.prfx}-val02_rare_v4.pdf")
# ggsave(val02.rare.p04, filename = fname.val02.p04, device=cairo_pdf, width=7) #width=5 