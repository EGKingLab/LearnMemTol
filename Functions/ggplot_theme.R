font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size),
  legend.text = element_text(size = font_size),
  legend.title=element_text(size = font_size),
  plot.title = element_text(size = font_size + 1))

library(scales)
show_col(hue_pal()(2))
LM_cols <- hue_pal()(2)
