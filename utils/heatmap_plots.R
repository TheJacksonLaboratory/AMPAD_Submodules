require(ggplot2)

heatmap <- function(mtx, row.annots=NA, col.annots=NA, value.name="Value") {
  
  rast.data <- as.data.frame(stack(mtx)) %>% dplyr::select(Rows=row, Cols=col, Value=4)
  if (!is.na(row.annots)) {
    rast.data <- merge(rast.data, row.annots, by.x="Rows", by.y="row.names")
    colnames(rast.data)[4] <- "Row.Annots"
  } else {
    rast.data$Row.Annots <- NA
  }
  if (!is.na(col.annots)) {
    rast.data <- merge(rast.data, col.annots, by.x="Cols", by.y="row.names")
    colnames(rast.data)[5] <- "Col.Annots"
  } else {
    rast.data$Col.Annots <- NA
  }
  
  p <- ggplot(data=rast.data) +
    geom_raster(aes(x=Cols, y=Rows, fill=Value)) +
    facet_grid(
      rows=vars(Row.Annots), #ifelse(any(!is.na(rast.data$Row.Annots)), vars(Row.Annots), vars(Row.Annots)), 
      cols=vars(Col.Annots), #ifelse(any(!is.na(rast.data$Col.Annots)), vars(Col.Annots), NA), 
      scales="free", space="free"
    ) +
    labs(fill=value.name) +
    scale_fill_gradient(low="#FFFFFF", high="#CC0000") +
    theme_minimal() +
    theme(
      panel.border=element_blank(),
      axis.line=element_line(size=0.25, color="black"),
      panel.grid=element_blank(),
      axis.title=element_blank(),
      axis.text.x=element_text(angle=270, hjust=0),
      strip.background=element_rect(color="black", size=0.5)
    )
  if (all(is.na(rast.data$Row.Annots))) {
    p <- p +
      theme(strip.text.y=element_blank(), strip.background.y=element_blank())
  }
  if (all(is.na(rast.data$Col.Annots))) {
    p <- p +
      theme(strip.text.x=element_blank(), strip.background.x=element_blank())
  }
  
  return(p)
}
