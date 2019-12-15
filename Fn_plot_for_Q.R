plor_for_Q <- function (Reference_gene,
                        Target_gene,
                        y_axis,
                        line1,
                        line2,
                        cols = c("2 KN" = "#E66101",
                                 "0.2 KN" = "#5E3C99"),
                        lines = c("Col-0" = "dotted",
                                  "nrt1.13" = "solid"),
                        cols_2KN_point = "#FDB863",
                        cols_0.2KN_point = "#B2ABD2",
                        cols_2KN_FT = "#E66101",
                        cols_0.2KN_FT = "#5E3C99",
                        save1)
{
  require (stringr)
  require (dplyr)
  require (ggplot2)
  require (magrittr)
  require (ggpubr)
  data_T <- na.omit (data_T)
  data_as <- data_T %>% group_by (line, KN, DAS) %>% summarise (AVG = mean (Concentration), STD = sd (Concentration)) # table of avg & std
  
  line_point_line1 <- paste (Target_gene, "-", Reference_gene, "_", line1, "_p.pdf")
  line_point_line2 <- paste (Target_gene, "-", Reference_gene, "_", line2, "_p.pdf")
  
  a <- as.numeric (((data_as[data_as$line == line1 & data_as$KN == "2 KN" & data_as$DAS == 40, "AVG"] - data_as[data_as$line == line1 & data_as$KN == "2 KN" & data_as$DAS == 32, "AVG"]) / (40 - 32) * (33 - 32)) + data_as[data_as$line == line1 & data_as$KN == "2 KN" & data_as$DAS == 32, "AVG"])
  b <- as.numeric (((data_as[data_as$line == line1 & data_as$KN == "0.2 KN" & data_as$DAS == 40, "AVG"] - data_as[data_as$line == line1 & data_as$KN == "0.2 KN" & data_as$DAS == 32, "AVG"]) / (40 - 32) * (37 - 32)) + data_as[data_as$line == line1 & data_as$KN == "0.2 KN" & data_as$DAS == 32, "AVG"])
  c <- as.numeric (((data_as[data_as$line == line2 & data_as$KN == "2 KN" & data_as$DAS == 40, "AVG"] - data_as[data_as$line == line2 & data_as$KN == "2 KN" & data_as$DAS == 32,"AVG"]) / (40 - 32) * (34.3 - 32)) + data_as[data_as$line == line2 & data_as$KN == "2 KN" & data_as$DAS == 32, "AVG"])
  d <- as.numeric (((data_as[data_as$line == line2 & data_as$KN == "0.2 KN" & data_as$DAS == 52, "AVG"] - data_as[data_as$line == line2 & data_as$KN == "0.2 KN" & data_as$DAS == 40, "AVG"]) / (52 - 40) * (47 - 40)) + data_as[data_as$line == line2 & data_as$KN == "0.2 KN" & data_as$DAS == 40, "AVG"])
  
  symnum.args <- list(cutpoints = c(0, 0.05, 1), symbols = c("*", "ns"))
  
  # Lineplot + scatterplot of line1
  # windows()
  ggline (data = data_T [data_T$line == line1, ],
          x = "DAS", y = "Concentration",
          add = "mean_sd",
          color = "KN") +
    
    labs (x = "DAS", y = y_axis,
          color = c("Nitrate concentration", "t.test"),
          shape = "Flowering time",
          title = line1) +
    
    scale_color_manual (values = cols,
                        breaks = c("2 KN", "0.2 KN")) +
    scale_linetype_manual (labels = c(line1),
                           values = lines) +
    scale_shape_manual (values = 4) +
    
    stat_compare_means (data = data_T [data_T$line == line1, ],
                        aes (group = KN),
                        method = "t.test",
                        label = "p.signif",
                        hide.ns = T,
                        color = "#1B9E77",
                        label.y = max (data_T$Concentration) * 1.05,
                        size = 12,
                        symnum.args = symnum.args) +
    
    geom_point (data = data_T [data_T$line == line1 & data_T$KN == "2 KN", ],
                aes (x = DAS, y = Concentration),
                size = 3,
                shape = 19,
                color = cols_2KN_point) +
    geom_point (data = data_T [data_T$line == line1 & data_T$KN == "0.2 KN", ],
                aes (x = DAS, y = Concentration),
                size = 3,
                shape = 19,
                color = cols_0.2KN_point) +
    
    geom_line (data = data_as [data_as$line == line1 & data_as$KN == "2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "2 KN"),
               size = 1.2) +
    geom_point (data = data_as [data_as$line == line1 & data_as$KN == "2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "2 KN"),
                size = 3,
                shape = 19) +
    geom_line (data = data_as [data_as$line == line1 & data_as$KN == "0.2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "0.2 KN"),
               size = 1.2) +
    geom_point (data = data_as [data_as$line == line1 & data_as$KN == "0.2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "0.2 KN"),
                size = 3,
                shape = 19) +
    
    geom_point (x = 33,
                y = a,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_2KN_FT,
                size = 1.5,
                stroke = 3) +
    geom_point (x = 37, y = b,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_0.2KN_FT,
                size = 1.5,
                stroke = 3) +
    
    theme (panel.background = element_rect (fill = "transparent"),
           legend.key = element_rect (color = "transparent",
                                      fill = "transparent"),
           legend.title = element_text (size = 24),
           legend.text = element_text (size = 20),
           legend.text.align = 0,
           legend.position = "right",
           legend.box.spacing = unit (1, "cm"),
           plot.title = element_text (size = 36,
                                      hjust = 0.5),
           axis.ticks.length.x = unit (0.3, "lines"),
           axis.line = element_line (color = "black",
                                     size = 1,
                                     linetype = "solid"),
           axis.text = element_text (size = 24,
                                     color = "black"),
           axis.title = element_text (size = 27),
           axis.title.x = element_text (vjust = -0.3),
           axis.title.y = element_text (vjust = 3),
           plot.margin = unit(c(0.8, 0.5, 0.8, 1), "cm"),
           panel.grid.major.y = element_line (color = "grey60")) +
    
    scale_x_continuous (breaks = c(6, 8, 10, 12, 14, 18, 21, 25, 28, 32, 40, 52)) +
    scale_y_continuous (limits = c(min (data_T$Concentration) * 0.9 , max (data_T$Concentration) * 1.1)) +
    
    ggsave (filename = line_point_line1,
            width = 15,
            height = 8,
            path = save1)
  
  # Lineplot + scatterplot of nrt1.13
  ggline (data = data_T [data_T$line == line2, ],
          x = "DAS", y = "Concentration",
          add = "mean_sd",
          color = "KN") +
    
    labs (x = "DAS", y = y_axis,
          color = c("Nitrate concentration", "t.test"),
          shape = "Flowering time",
          title = line2) +
    
    scale_color_manual (values = cols,
                        breaks = c("2 KN", "0.2 KN")) +
    scale_linetype_manual ("Line",
                           labels = c(line1, expression(italic(line2))),
                           values = lines) +
    scale_shape_manual (values = 4) +
    
    stat_compare_means (data = data_T [data_T$KN == "2 KN", ],
                        aes (group = line),
                        method = "t.test",
                        label = "p.signif",
                        hide.ns = T,
                        color = "#E66101",
                        label.y = max (data_T$Concentration) * 1.05,
                        size = 12,
                        symnum.args = symnum.args) +
    stat_compare_means (data = data_T [data_T$KN == "0.2 KN", ],
                        aes (group = line),
                        method = "t.test",
                        label = "p.signif",
                        hide.ns = T,
                        color = "#5E3C99",
                        label.y = max (data_T$Concentration) * 1.1,
                        size = 12,
                        symnum.args = symnum.args) +
    stat_compare_means (data = data_T [data_T$line == line1, ],
                        aes (group = KN),
                        method = "t.test",
                        label = "p.signif",
                        hide.ns = T,
                        color = "#1B9E77",
                        label.y = max (data_T$Concentration) * 1.15,
                        size = 12,
                        symnum.args = symnum.args) +
    stat_compare_means (data = data_T [data_T$line == line2, ],
                        aes (group = KN),
                        method = "t.test",
                        label = "p.signif",
                        hide.ns = T,
                        color = "#E7298A",
                        label.y = max (data_T$Concentration) * 1.2,
                        size = 12,
                        symnum.args = symnum.args) +
    
    geom_line (data = data_as [data_as$line == line1 & data_as$KN == "2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "2 KN",
                    linetype = line1),
               size = 1.2) +
    geom_point (data = data_as [data_as$line == line1 & data_as$KN == "2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "2 KN"),
                size = 3,
                shape = 21) +
    geom_line (data = data_as [data_as$line == line1 & data_as$KN == "0.2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "0.2 KN",
                    linetype = line1),
               size = 1.2) +
    geom_point (data = data_as [data_as$line == line1 & data_as$KN == "0.2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "0.2 KN"),
                size = 3,
                shape = 21) +
    
    geom_point (data = data_T [data_T$line == line2 & data_T$KN == "2 KN", ],
                aes (x = DAS, y = Concentration),
                size = 3,
                shape = 19,
                color = cols_2KN_point) +
    geom_point (data = data_T [data_T$line == line2 & data_T$KN == "0.2 KN", ],
                aes (x = DAS, y = Concentration),
                size = 3,
                shape = 19,
                color = cols_0.2KN_point) +
    
    geom_point (data = data_as [data_as$line == line2 & data_as$KN == "2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "2 KN"),
                size = 3,
                shape = 19) +
    geom_point (data = data_as [data_as$line == line2 & data_as$KN == "0.2 KN", ],
                aes (x = DAS, y = AVG,
                     color = "0.2 KN"),
                size = 3,
                shape = 19) +
    
    geom_line (data = data_as [data_as$line == line2 & data_as$KN == "2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "2 KN",
                    linetype = line2),
               size = 1.2) +
    geom_line (data = data_as [data_as$line == line2 & data_as$KN == "0.2 KN", ],
               aes (x = DAS, y = AVG,
                    color = "0.2 KN",
                    linetype = line2),
               size = 1.2) +
    
    geom_point (x = 33, y = a,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_2KN_FT,
                size = 1.5,
                stroke = 3) +
    geom_point (x = 37, y = b,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_0.2KN_FT,
                size = 1.5,
                stroke = 3) +
    geom_point (x = 34.3, y = c,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_2KN_FT,
                size = 2,
                stroke = 3) +
    geom_point (x = 47, y = d,
                aes (x = DAS, y = AVG,
                     shape = "FT"),
                color = cols_0.2KN_FT,
                size = 2,
                stroke = 3) +
    
    theme (panel.background = element_rect (fill = "transparent"),
           legend.key = element_rect (color = "transparent",
                                      fill = "transparent"),
           legend.title = element_text (size = 24),
           legend.text = element_text (size = 20),
           legend.text.align = 0,
           legend.position = "right",
           legend.box.spacing = unit (1, "cm"),
           plot.title = element_text (size = 36,
                                      face = "italic",
                                      hjust = 0.5),
           axis.ticks.length.x = unit (0.3, "lines"),
           axis.line = element_line (color = "black",
                                     size = 1,
                                     linetype = "solid"),
           axis.text = element_text (size = 24,
                                     color = "black"),
           axis.title = element_text (size = 27),
           axis.title.x = element_text (vjust = -0.3),
           axis.title.y = element_text (vjust = 3),
           plot.margin = unit(c(0.8, 0.5, 0.8, 1), "cm"),
           panel.grid.major.y = element_line (color = "grey60")) +
    
    scale_x_continuous (breaks = c(6, 8, 10, 12, 14, 18, 21, 25, 28, 32, 40, 52)) +
    scale_y_continuous (limits = c(min (data_T$Concentration) * 0.9 , max (data_T$Concentration) * 1.3)) +
    
    ggsave (filename = line_point_line2,
            width = 15,
            height = 8,
            path = save1)
}
