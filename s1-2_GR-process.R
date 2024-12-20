
# Calculate GRs using GrowthRateR
#EXP20220601.gr <- growthrateR(EXP20220601.raw, EXP20220601.platemap, time.max=5, plate.clean=TRUE)
EXP20220602.gr <- growthrateR(EXP20220602.raw, EXP20220602.platemap, time.max=5, plate.clean=TRUE)
EXP20220605.gr <- growthrateR(EXP20220605.raw, EXP20220605.platemap, time.max=5, plate.clean=TRUE)
#EXP20220625A.gr <- growthrateR(EXP20220625A.raw, EXP20220625A.platemap, time.max=5, plate.clean=TRUE)
EXP20220626.gr <- growthrateR(EXP20220626.raw, EXP20220626.platemap, time.max=5, plate.clean=TRUE)
#EXP20220627.gr <- growthrateR(EXP20220627.raw, EXP20220627.platemap, time.max=5, plate.clean=TRUE)

# Get list of exps want to keep
gr_list <- list(#EXP20220601.gr,
  EXP20220602.gr, 
  EXP20220605.gr,
  #EXP20220625A.gr,
  EXP20220626.gr
  #EXP20220627.gr
)

# Merge together - this will have the individual replicates
gr.merged <- do.call(rbind, gr_list)
gr.merged <- gr.merged %>%
  mutate(code = as.character(code)) %>%
  mutate(test = paste(code, exp)) %>%
  na.omit()

# Get summary stats
gr.merged.stats <- gr.merged %>%
  group_by(code, gene, exp, test, IPTG) %>%
  summarise(N=length(slope), Average=mean(slope), StDev=sd(slope), AvgTime=mean(time)) %>%
  ungroup() %>%
  na.omit()

# Make list of duplicate tests to remove
filter_list <- c("351 EXP20220602",
                 "340 EXP20220626",
                 #"342 EXP20220626",
                 "344 EXP20220602",
                 "344 EXP20220605",
                 "345 EXP20220626",
                 "346 EXP20220602",
                 "351 EXP20220602",
                 "352 EXP20220602",
                 "353 EXP20220602",
                 #"359 EXP20220605",
                 #"360 EXP20220605",
                 #"361 EXP20220605",
                 "363 EXP20220605",
                 "363 EXP20220626"
)

# Remove 'em from both dfs
gr.merged.fil <- gr.merged %>% filter(!test %in% filter_list)
gr.merged.stats.fil <- gr.merged.stats %>% filter(!test %in% filter_list)


# Export file for Prism analyses
gr.merged.noIPTG <- filter(gr.merged.fil, IPTG != "N")
gr.fname <- stringr::str_interp("${out.path}/${fln.prfx}-rep_raw.csv")
#write.csv(file=gr.fname, gr.merged.fil, quote=F, row.names=FALSE)
write.csv(file=gr.fname, gr.merged.noIPTG, quote=F, row.names=FALSE)

# Export so don't have to calculate each time
gr.merged.fil.noIPTG <- filter(gr.merged.stats.fil, IPTG != "N")
gr.noIPTG.fname <- stringr::str_interp("${out.path}/${fln.prfx}-gr_merged.csv")
write.csv(file=gr.noIPTG.fname, gr.merged.fil.noIPTG, quote=F, row.names=FALSE)
