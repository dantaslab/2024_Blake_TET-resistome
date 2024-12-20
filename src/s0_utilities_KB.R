### Set Colors -----------------------------------------------------------------

log2.pal <- c("0" = "#FFFFFF",
              "1" = "#EEDEDE",
              "2" = "#DDBDBD",
              "3" = "#CD9C9C",
              "4" = "#BC7C7C",            
              "5" = "#AC5B5B",
              "6" = "#9B3A3A" #,
              #"7" = "#8B1A1A"
)

log2.color.count <- 8
log2.colors <- c("white", "firebrick4")
log2.colors.pal <- colorRampPalette(log2.colors)(8)

mech.pal2 <- list(mech = c(
            "DES1" = "#7570B3",
            "DES2" = "#E35231",
            "EFF" = "#1B9E77",
            "none" = "#4D4D4D",
            "RPP" = "#FBCB6A"))

mech.pal <- c(
            "DES1" = "#7570B3",
            "DES2" = "#E35231",
            "EFF" = "#1B9E77",
            "none" = "gray50",
            "RPP" = "#FBCB6A")

gene.pal <- c(#DES1 - Blues
              "tetX2" = "#2f2c53",
              "tetX3" = "#423e75", 
              "tetX6" = "#555096", 
              "tetX7" = "#6e69af", 
              "tetX8" = "#8f8ac1",
              "tetX12" = "#afacd3", 

              #EFF - Greens
              "tetA" = "#136d52",
              "tetB" = "#1a9873",
              "tetE" = "#21c494",
              "tetG" = "#3bdead",
              "tetL" = "#67e5bf",
              "tet39" = "#92ecd2",
              
              #DES2 - Reds
              "tet47" = "#70210f",
              "tet50" = "#9d2f15",
              "tet51" = "#ca3c1b",
              "tet54" = "#e45535",
              "tet55" = "#ea7b62",
              "tet56" = "#f0a18f",
              
              #RPP - Yellows - UPDATE!!!
              "tetM" = "#7c4d04",
              "tetO" = "#ad6c05",
              "tetS" = "#df8b07",
              "tetW" = "#f8a520",
              "tet32" = "#fab952",
              "tet36" = "#fbcd83",

              "empty" = "gray50"
)

### General orderings ----------------------------------------------------------

DES1.list <- c("tetX2", "tetX3", "tetX6", "tetX7", "tetX8", "tetX12")
EFF.list <- c("tetA", "tetB", "tetE", "tetG", "tetL", "tet39")
RPP.list <- c("tetM", "tetO", "tetS", "tetW", "tet32", "tet36")
DES2.list <- c("tet47", "tet50", "tet51", "tet54", "tet55", "tet56")

TET.order <- c("none-none-0", "TET-0.5-20", "TET-1-20", "TET-4-20")
DOX.order <- c("none-none-0", "DOX-0.5-20", "DOX-1-20", "DOX-4-20")
MIN.order <- c("none-none-0", "MIN-0.5-20", "MIN-1-20", "MIN-4-20")
TIG.order <- c("none-none-0", "TIG-0.25-20", "TIG-0.5-20", "TIG-2-48")
none.order <- c("none-none-0", "none-none-20")
ATC.order <- c("none-none-0", "ATC-2-24", "ATC-4-24", "ATC-8-24")

DES2.TET.order <- c("none-none-0", "TET-0.5-20", "TET-1-20", "TET-4-48")
DES2.DOX.order <- c("none-none-0", "DOX-0.5-20", "DOX-1-20", "DOX-4-48")
DES2.MIN.order <- c("none-none-0", "MIN-0.5-20", "MIN-1-20", "MIN-4-48")
DES1.TIG.order <- c("none-none-0", "TIG-0.25-20", "TIG-0.5-48", "TIG-2-48")
RPP.MIN.order <- c("none-none-0", "MIN-0.5-20", "MIN-1-20", "MIN-4-48")
EFF.MIN.order <- c("none-none-0", "MIN-0.5-20", "MIN-1-20", "MIN-4-48")
minDES2.TIG.order <- c("none-none-0", "TIG-0.25-20", "TIG-0.5-20", "TIG-2-48") #0.5-48

plot.name <- c("none-none-0" = "input mix",
               "TET-0.5-20" = "TET 0.5x",
               "TET-1-20" = "TET 1x",
               "TET-4-20" = "TET 4x",
               "DOX-0.5-20" = "DOX 0.5x",
               "DOX-1-20" = "DOX 1x",
               "DOX-4-20" = "DOX 4x",
               "MIN-0.5-20" = "MIN 0.5x",
               "MIN-1-20" = "MIN 1x",
               "MIN-4-20" = "MIN 4x",
               "TIG-0.25-20" = "TIG 0.5x",
               "TIG-0.5-20" = "TIG 1x",
               "TIG-2-48" = "TIG 4x"
               )

### Helper functions -----------------------------------------------------------

export_plot <- function(plot, fname, w, h, fmt=c("pdf", "png")) {
  
  fname.pdf <- stringr::str_interp("${out.path}/${fln.prfx}-${fname}.pdf")
  fname.png <- stringr::str_interp("${out.path}/${fln.prfx}-${fname}.png")
  
  if (missing(fmt)){
    fmt <- "pdf"
  }
  if (fmt == "pdf") {
    save <- ggsave(plot, filename=fname.pdf, device=cairo_pdf, width=w, height=h)
  }
  if (fmt == "png") {
    save <- ggsave(plot, filename=fname.png, device=png, width=w, height=h)
  }
  
  return(save)
}


### Custom themes --------------------------------------------------------------

theme_pub <- function() {

  color.background = "white"
  color.grid.major = "grey30"
  color.grid.minor = "grey60"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = "black"
  
  # Begin construction of chart
    theme_cowplot() +
    ggplot2::theme( 
        text = element_text(family = "Arial", size=6),
        line = element_line(size=0.3)
    ) +
      
    # Format the legend, hide by default
    ggplot2::theme(  
        legend.position    = "bottom",
        #legend.position    = "none",
        legend.margin      = margin(4, 4, 4, 4),
        legend.key.size    = unit(8, "pt"),
        legend.box.spacing = unit(4, "pt"),
        legend.text        = element_text(size=6, color = color.axis.title),
        legend.title       = element_text(size=8, face="bold")
    ) +
  
    # Set axis labels, text and format tick marks
    ggplot2::theme(  
           axis.text          = element_text(size=6),
           axis.text.x        = element_text(size = 6, angle = 0, color=color.axis.text, margin = margin(1, 0, 0, 0)),
           axis.text.y        = element_text(size = 6, angle = 0, color=color.axis.text, margin = margin(0, 1, 0, 0)),
           axis.text.x.top    = element_text(margin = margin(0, 0, 1, 0)),
           axis.text.y.right  = element_text(margin = margin(0, 0, 0, 1)),
           axis.ticks         = element_line(size=0.3),
           axis.ticks.length  = unit(2, "pt"),

           axis.line          = element_line(size=0.3),
           axis.line.x        = element_line(size=0.3),
           axis.line.y        = element_line(size=0.3),
            
           axis.title         = element_text(size=10, face="bold", color = color.axis.title),
           axis.title.x.top   = element_text(margin = margin(0, 0, 2, 0)),
           axis.title.y.right = element_text(margin = margin(0, 0, 0, 2)),
           axis.title.x       = element_text(margin = margin(1, 0, 0.5, 0), vjust=0),
           axis.title.y       = element_text(margin = margin(0, 1, 0, 0.5), vjust=1.25)
    ) +

    # Format facet heading
    ggplot2::theme( 
            strip.text   = element_text(size=10, face="bold"),
            strip.text.x = element_text(margin = margin(3, 0, 1.5, 0)), 
            strip.text.y = element_text(margin = margin(0, 3, 0, 3)),
            # strip.background = element_rect(fill = "#efefef", color = "#2a3132"),
            strip.background   = element_blank()
    ) +
        
    ## Format the grid
    ggplot2::theme(
    #        panel.grid.major=element_line(color=color.grid.major,size=.50)) +
    #        panel.grid.minor=element_line(color=color.grid.minor,size=.50)) +
             panel.spacing      = unit(3, "pt")
    ) +
    
    # Format plot title and margin  
    ggplot2::theme(  
            plot.title = element_text(size=10, color=color.title),
            # plot.title = element_text(size=12, color=color.title, family = "Futura-Book"),
            # plot.margin = unit(c(0.35, 0.3, 0.35, 0.35), "cm"),
            plot.margin = margin(3, 3, 3, 3)
            
    )
}

theme_plot <- function() {
  
  color.background = "white"
  color.grid.major = "grey30"
  color.grid.minor = "grey60"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = "black"
  
  # Begin construction of chart
  theme_cowplot() +
    ggplot2::theme( 
      text = element_text(family = "Arial", size=6),
      line = element_line(size=0.3)
    ) +
    
    # Format the legend, hide by default
    ggplot2::theme(  
      legend.position    = "none",
      legend.margin      = margin(4, 4, 4, 4),
      legend.key.size    = unit(8, "pt"),
      legend.box.spacing = unit(4, "pt"),
      legend.text        = element_text(size=6, color = color.axis.title),
      legend.title       = element_text(size=8, face="bold")
    ) +
    
    # Set axis labels, text and format tick marks
    ggplot2::theme(  
      axis.text          = element_text(size=6),
      axis.text.x       = element_text(size=6, angle = 45, vjust=1, hjust=1), #, hjust=1, vjust=0.5
      axis.text.y        = element_text(size = 6, angle = 0, color=color.axis.text, margin = margin(0, 1, 0, 0)),
      axis.text.x.top    = element_text(margin = margin(0, 0, 1, 0)),
      axis.text.y.right  = element_text(margin = margin(0, 0, 0, 1)),
      axis.ticks         = element_line(size=0.75),
      axis.ticks.length  = unit(2, "pt"),
      
      #axis.line          = element_line(size=0.3),
      axis.line.x        = element_blank(),
      axis.line.y        = element_blank(),
      panel.border       = element_rect(colour = "black", fill=NA, size=1),
      
      axis.title         = element_text(size=8, face="bold", color = color.axis.title),
      axis.title.x.top   = element_text(margin = margin(0, 0, 2, 0)),
      axis.title.y.right = element_text(margin = margin(0, 0, 0, 2)),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(margin = margin(0, 1, 0, 0.5), vjust=1.25)
    ) +
    
    # Format facet heading
    ggplot2::theme( 
      strip.text   = element_text(size=8, face="bold"),
      strip.text.x = element_text(margin = margin(3, 0, 1.5, 0)), 
      strip.text.y = element_text(margin = margin(0, 3, 0, 3)),
      # strip.background = element_rect(fill = "#efefef", color = "#2a3132"),
      strip.background   = element_blank()
    ) +
    
    ## Format the grid
    ggplot2::theme(
      #        panel.grid.major=element_line(color=color.grid.major,size=.50)) +
      #        panel.grid.minor=element_line(color=color.grid.minor,size=.50)) +
      panel.spacing      = unit(1, "pt") #3
    ) +
    
    # Format plot title and margin  
    ggplot2::theme(  
      plot.title = element_text(size=10, color=color.title),
      # plot.title = element_text(size=12, color=color.title, family = "Futura-Book"),
      # plot.margin = unit(c(0.35, 0.3, 0.35, 0.35), "cm"),
      plot.margin = margin(1, 1, 1, 1) #3
      
      
    )
}
