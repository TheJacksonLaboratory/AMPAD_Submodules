require(ggplot2)

basic_theme <- theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    axis.line=element_line(color="black")
  )

basic_theme_big <- theme_linedraw(base_size=18) +
  basic_theme
