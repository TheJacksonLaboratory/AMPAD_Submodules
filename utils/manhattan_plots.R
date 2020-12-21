# Manhattan Plots

manhattan.plot <- function(gwas.results, regions=data.frame(), gwsug=1e-05, gwsig=5e-08) {
  
  man.plot <- gwas.results %>% 
    dplyr::group_by(CHR) %>% 
    dplyr::summarise(chr_len=max(BP)) %>% 
    dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    dplyr::left_join(gwas.results, ., by=c("CHR"="CHR")) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum=BP+tot)
  
  axisdf = man.plot %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  p <- ggplot(man.plot, aes(x=BPcum, y=-log10(P))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    geom_hline(yintercept = -log10(gwsig), color="firebrick1", lty=2) +
    geom_hline(yintercept = -log10(gwsug), color="goldenrod1", lty=2) +
    scale_color_manual(values = rep(c("royalblue1", "skyblue"), 22)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(-log10(5e-02), 12), breaks=2:12 ) +
    labs(x="", y="-log10(p)") +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (nrow(regions) == 0) {
    return(p)
  }
  
  for (i in 1:nrow(regions)) {
    p <- p + geom_point(
      data=man.plot %>% dplyr::filter(CHR == regions[i,"CHR"]) %>%
        dplyr::filter(BP <= regions[i,"BP.End"] & BP >= regions[i,"BP.Start"]),
      mapping=aes(x=BPcum, y=-log10(P)),
      color=regions[i,"Color"],
      alpha=0.8, size=1.2
    )
  }
  
  return(p)
}
