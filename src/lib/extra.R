test <- data.frame(sample = paste0("sample", 1:100), 
    tumor = rep(c("breast", "lung", "bladder", "colorectal", "prostate"), 20), 
    allele = rep(c("exon1", "exon2", "exon3", "exon4", "exon5", "exon6", "exon7", "exon8", "exon9", "exon10"),10), 
    values = rnorm(100))
test <- test[order(test$tumor),]
test$sample <- factor(test$sample, levels = test$sample)

h1 <- ggplot(test)+
    geom_bar(mapping = aes(x = sample, y = values, fill = allele), 
        stat = "identity", 
        position = "dodge")+
    guides(fill = guide_legend(ncol = 2))+
    theme(axis.text.x.bottom = element_blank(), 
        axis.ticks = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank())+
    facet_grid(.~tumor, scales = "free_x")
h2 <- ggplot(test)+
    geom_bar(mapping = aes(x = sample, y = 1, fill = tumor), 
        stat = "identity", 
        width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))+
    facet_grid(.~tumor, scales = "free_x")

legend <- cowplot::plot_grid(get_legend(h2), get_legend(h1), ncol = 1)
h1 <- h1 + theme(legend.position = "none")
h2 <- h2 + theme(legend.position = "none")
plot <- cowplot::plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 15))
cowplot::plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.5))

allele_col = structure(1:10, names = paste0("exon", 1:10))
ha_list <- HeatmapAnnotation(bar = anno_barplot(test$values,
    gp = gpar(fill = allele_col[test$allele]),
    height = unit(6, "cm"))) %v%
    HeatmapAnnotation(tumor = anno_simple(test$tumor))

draw(ha_list, heatmap_legend_list = Legend(title = "Allele", labels = names(allele_col), legend_gp = gpar(fill = allele_col)))