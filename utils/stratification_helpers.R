# Goal
# * Develop functions to plot module stratification
# * Develop a cluster decision making algorithm

# Load libraries
library(plotrix)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# Plots a matrix-like object in a colored heatmap
# Matrix object must have colnames and rownames defined
# Matrix object must contain only double values
plot.strat <- function(mtx, main="Stratification", xlab="Module", ylab="Cluster", text=T) {
  
  plot.df <- mtx %>%
    as.data.frame()
  plot.df$ID <- rownames(plot.df)
  plot.df <- melt(plot.df, id.vars="ID")
  plot.df$ID <- factor(plot.df$ID, levels=unique(plot.df$ID))
  
  p <- ggplot(data=plot.df) + 
    geom_raster(aes(variable, ID, fill=value)) +
    geom_text(aes(variable, ID, label=round(value, 2))) +
    geom_hline(yintercept=(0:nrow(mtx))+0.5, color="#444444", lty=1) +
    geom_vline(xintercept=(0:ncol(mtx))+0.5, color="#444444", lty=1) +
    scale_fill_distiller(palette="RdYlBu") +
    scale_x_discrete() +
    scale_y_discrete() +
    guides(fill=guide_colorbar(title="Subtype")) +
    theme_linedraw(base_size=18) + 
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(main) +
    theme(axis.text.x=element_text(angle=90, hjust=1))
  
  # For poster only
  if (!text) {
    p <- ggplot(data=plot.df) + 
      geom_raster(aes(variable, ID, fill=value)) +
      geom_hline(yintercept=(0:nrow(mtx))+0.5, color="black", lty=1) +
      geom_vline(xintercept=(0:ncol(mtx))+0.5, color="black", lty=1) +
      scale_fill_distiller(palette="RdYlBu") +
      guides(fill=guide_colorbar(title="Scale")) +
      scale_x_discrete() +
      scale_y_discrete() +
      theme_linedraw(base_size=18) +
      xlab(xlab) +
      ylab(ylab) + 
      ggtitle(main) +
      theme(axis.text.x=element_text(angle=60, hjust=1), axis.title.x=element_blank())
  }

  return(p)
}

# Plots a matrix-like object in a colored heatmap
# Matrix object must have colnames and rownames defined
# Matrix object must contain only double values
# Used for large number of values
plot.strat.large <- function(mtx, main="Stratification", xlab="Module", ylab="Cluster") {
  
  plot.df <- mtx %>%
    as.data.frame()
  plot.df$ID <- rownames(plot.df)
  plot.df <- melt(plot.df, id.vars="ID")
  plot.df$ID <- factor(plot.df$ID, levels=unique(plot.df$ID))
  
  p <- ggplot(data=plot.df) +
    geom_raster(aes(variable, ID, fill=value)) +
    scale_fill_gradient2(low="#3300FF", high="#FF0033") +
    guides(fill=guide_colorbar(title="Scale")) +
    scale_x_discrete() +
    scale_y_discrete() +
    theme_linedraw(base_size=18) +
    xlab(xlab) + 
    ylab(ylab) +
    ggtitle(main) +
    theme(axis.text=element_blank(), axis.ticks=element_blank())
  
  return(p)
}

# Vector must has the same cardinality as the number of patients
# and must be named with the patient IDs from the cleaned data.
# This is specifically for 12_stratification.R
trait.dist <- function(trait, name, categorical) {
  
  # Get method-specific trait values
  trait.control <- trait[eigengenes.controls$Patient]
  
  # Get a list of distributions of the trait for each cluster
  dists <- sapply(
    1:length(levels(pca.cluster.case$Best.partition)), 
    function(x) trait[names(pca.cluster.case$Best.partition[pca.cluster.case$Best.partition == x])]
  )
  
  # Organize into a data frame
  dists.tbl <- data.frame(
    Trait=trait,
    Population="Total"
  )
  dists.tbl <- rbind(dists.tbl, data.frame(
    Trait=trait.control,
    Population="Control"
  ))
  dists.tbl <- rbind(
    dists.tbl,
    do.call("rbind", lapply(seq_along(dists), function(x, i) data.frame(Trait=x[[i]], Population=LETTERS[i]), x=dists))
  )
  
  dists.tbl <- na.omit(dists.tbl)
  
  # Generate a jitter plot to view effect
  p <- ggplot() +
    geom_jitter(data=dists.tbl, aes(x=Population, y=Trait, color=Population), size=I(2), alpha=0.75, width=0.2, height=0) +
    ggtitle(paste0(name, " Distribution Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Trait Value") +
    theme_linedraw(base_size=18)
  
  q <- ggplot() +
    geom_bar(data=dists.tbl, aes(x=Population, fill=Trait), position=position_fill(), color="black") +
    ggtitle(paste0(name, " Proportions Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Proportion") +
    theme_linedraw(base_size=18)
  
  r <- ggplot() +
    geom_boxplot(data=dists.tbl, aes(x=Population, y=Trait, fill=Population)) +
    xlab("Population") + ylab("Trait Value") + 
    guides(fill=F) +
    scale_fill_brewer(palette="Set3") +
    theme_linedraw(base_size=18) +
    theme(
      plot.title=element_blank(),
      plot.subtitle=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.border=element_blank(),
      axis.line=element_line(color="black")
    )  
  
  # Add means to continuous variables
  if (!categorical) {
    means <- as.data.frame(dists.tbl %>% group_by(Population) %>% summarise(mean = mean(Trait)))
    rownames(means) <- means$Population
    means <- means[dists.tbl$Population,]
    p <- p + geom_point(data=means, aes(x=Population, y=mean), pch=4, size=I(4))
  }
  
  # Run Single-Factor ANOVA
  dists.tbl$Trait = as.numeric(dists.tbl$Trait)
  res.aov <- aov(Trait ~ Population, data=subset(dists.tbl, Population != "Control"))
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  
  # If categorical, run a Chi-Square test
  # If not, run a T-Test
  if (categorical) {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      c.tbl <- table(dists.tbl$Trait, dists.tbl$Population)[,c("Total", LETTERS[i])]
      p.vals <- c(p.vals, chisq.test(c.tbl)$p.value)
    }
  } else {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      p.vals <- c(p.vals, t.test(subset(dists.tbl, Population=="Total")$Trait, subset(dists.tbl, Population==LETTERS[i])$Trait)$p.value)
    }
  }
  
  p <- p + geom_vline(xintercept=(2+1:length(levels(pca.cluster.case$Best.partition)))[p.vals < 0.05], lty=2, size=1)
  r <- r + geom_text(aes(x=levels(dists.tbl$Population), y=sapply(split(dists.tbl, dists.tbl$Population), function(x) mean(x$Trait)), label=c("", "", ifelse(p.vals < 0.05, "*", ""))), size=8, vjust=-0.05)
  
  # Create bar graph
  if (categorical) {
    q <- q + geom_text(
      aes(x=1+1:length(levels(pca.cluster.case$Best.partition)), y=1, label=ifelse(p.vals < 0.05, "*", "")), 
      vjust=0.2, color="black", position=position_dodge(0.9), size=7.5
    )
    return(list(dist.plot=p, prop.plot=q))
  }
  
  return(list(dist.plot=p, box.plot=r))
}

# Vector must has the same cardinality as the number of patients
# and must be named with the patient IDs from the cleaned data
# This is specifically for 12a_stratification.R
trait.dist.a <- function(trait, name, categorical=T) {
  
  # Get method-specific trait values
  trait.control <- trait[rownames(eigengenes.controls)]
  trait.other <- trait[rownames(eigengenes.other)]
  
  # Get a list of distributions of the trait for each cluster
  dists <- sapply(
    1:length(levels(pca.cluster.case$Best.partition)), 
    function(x) trait[names(pca.cluster.case$Best.partition[pca.cluster.case$Best.partition == x])]
  )
  
  # Organize into a data frame
  dists.tbl <- data.frame(
    Trait=trait,
    Population="Total"
  )
  dists.tbl <- rbind(dists.tbl, data.frame(
    Trait=trait.control,
    Population="Control"
  ))
  dists.tbl <- rbind(dists.tbl, data.frame(
    Trait=trait.other,
    Population="Other"
  ))
  dists.tbl <- rbind(
    dists.tbl,
    do.call("rbind", lapply(seq_along(dists), function(x, i) data.frame(Trait=x[[i]], Population=LETTERS[i]), x=dists))
  )
  
  dists.tbl <- na.omit(dists.tbl)
  
  # Generate a jitter plot to view effect
  p <- ggplot() +
    geom_jitter(data=dists.tbl, aes(x=Population, y=Trait, color=Population), size=I(2), alpha=0.75, width=0.2, height=0) +
    ggtitle(paste0(name, " Distribution Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Trait Value") +
    theme_linedraw(base_size=18)
  
  q <- ggplot() +
    geom_bar(data=dists.tbl, aes(x=Population, fill=Trait), position=position_fill(), color="black") +
    ggtitle(paste0(name, " Proportions Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Proportion") +
    theme_linedraw(base_size=18)
  
  # Add means to continuous variables
  if (!categorical) {
    means <- as.data.frame(dists.tbl %>% group_by(Population) %>% summarise(mean = mean(Trait)))
    rownames(means) <- means$Population
    means <- means[dists.tbl$Population,]
    p <- p + geom_point(data=means, aes(x=Population, y=mean), pch=4, size=I(4))
  }
  
  # Run Single-Factor ANOVA
  dists.tbl$Trait = as.numeric(dists.tbl$Trait)
  res.aov <- aov(Trait ~ Population, data=subset(dists.tbl, Population != "Other" & Population != "Control"))
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  
  # If categorical, run a Chi-Square test
  # If not, run a T-Test
  if (categorical) {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      c.tbl <- table(dists.tbl$Trait, dists.tbl$Population)[,c("Total", LETTERS[i])]
      p.vals <- c(p.vals, chisq.test(c.tbl)$p.value)
    }
  } else {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      p.vals <- c(p.vals, t.test(subset(dists.tbl, Population=="Total")$Trait, subset(dists.tbl, Population==LETTERS[i])$Trait)$p.value)
    }
  }
  
  p <- p + geom_vline(xintercept=(3+1:length(levels(pca.cluster.case$Best.partition)))[p.vals < 0.05], lty=2, size=1)
  
  # Create bar graph
  if (categorical) {
    q <- q + geom_text(
      aes(x=1+1:length(levels(pca.cluster.case$Best.partition)), y=1, label=ifelse(p.vals < 0.05, "*", "")), 
      vjust=0.2, color="black", position=position_dodge(0.9), size=7.5
    )
    return(list(dist.plot=p, prop.plot=q))
  }
  
  return(p)
}

# Vector must has the same cardinality as the number of patients
# and must be named with the patient IDs from the cleaned data
# This is specifically for 12b_stratification.R
trait.dist.b <- function(trait, name, categorical=T, population=LETTERS) {
  
  # Get a list of distributions of the trait for each cluster
  dists <- sapply(
    1:length(levels(pca.cluster.case$Best.partition)), 
    function(x) trait[names(pca.cluster.case$Best.partition[pca.cluster.case$Best.partition == x])]
  )
  
  # Organize into a data frame
  dists.tbl <- data.frame(
    Trait=trait,
    Population="Total"
  )
  dists.tbl <- rbind(
    dists.tbl,
    do.call("rbind", lapply(seq_along(dists), function(x, i) data.frame(Trait=x[[i]], Population=population[i]), x=dists))
  )
  
  dists.tbl <- na.omit(dists.tbl)
  
  # Generate a jitter plot to view effect
  p <- ggplot() +
    geom_jitter(data=dists.tbl, aes(x=Population, y=Trait, color=Population), size=I(2), alpha=0.75, width=0.2, height=0) +
    ggtitle(paste0(name, " Distribution Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Trait Value") +
    theme_linedraw(base_size=18)
  
  q <- ggplot() +
    geom_bar(data=dists.tbl, aes(x=Population, fill=Trait), position=position_fill(), color="black") +
    ggtitle(paste0(name, " Proportions Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Proportion") +
    theme_linedraw(base_size=18)
  
  # Add means to continuous variables
  if (!categorical) {
    means <- as.data.frame(dists.tbl %>% group_by(Population) %>% summarise(mean = mean(Trait)))
    rownames(means) <- means$Population
    means <- means[dists.tbl$Population,]
    p <- p + geom_point(data=means, aes(x=Population, y=mean), pch=4, size=I(4))
  }
  
  # Run Single-Factor ANOVA
  dists.tbl$Trait = as.numeric(dists.tbl$Trait)
  res.aov <- aov(Trait ~ Population, data=subset(dists.tbl, Population != "Control"))
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  
  # If categorical, run a Chi-Square test
  # If not, run a T-Test
  if (categorical) {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      c.tbl <- table(dists.tbl$Trait, dists.tbl$Population)[,c("Total", population[i])]
      p.vals <- c(p.vals, chisq.test(c.tbl)$p.value)
    }
  } else {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      p.vals <- c(p.vals, t.test(subset(dists.tbl, Population=="Total")$Trait, subset(dists.tbl, Population==population[i])$Trait)$p.value)
    }
  }
  
  p <- p + geom_vline(xintercept=(1+1:length(levels(pca.cluster.case$Best.partition)))[p.vals < 0.05], lty=2, size=1)
  
  # Create bar graph
  if (categorical) {
    q <- q + geom_text(
      aes(x=1+1:length(levels(pca.cluster.case$Best.partition)), y=1, label=ifelse(p.vals < 0.05, "*", "")), 
      vjust=0.2, color="black", position=position_dodge(0.9), size=7.5
    )
    return(list(dist.plot=p, prop.plot=q))
  }
  
  return(p)
}

# Vector must has the same cardinality as the number of patients
# and must be named with the patient IDs from the cleaned data
# This is specifically for 12c_stratification.R
trait.dist.c <- function(trait, name, categorical=T) {
  
  # Get a list of distributions of the trait for each cluster
  dists <- sapply(
    1:length(levels(pca.cluster.case$Best.partition)), 
    function(x) trait[names(pca.cluster.case$Best.partition[pca.cluster.case$Best.partition == x])]
  )
  
  # Organize into a data frame
  dists.tbl <- data.frame(
    Trait=trait,
    Population="Total"
  )
  dists.tbl <- rbind(
    dists.tbl,
    do.call("rbind", lapply(seq_along(dists), function(x, i) data.frame(Trait=x[[i]], Population=LETTERS[i]), x=dists))
  )
  
  dists.tbl <- na.omit(dists.tbl)
  
  # Generate a jitter plot to view effect
  p <- ggplot() +
    geom_jitter(data=dists.tbl, aes(x=Population, y=Trait, color=Population), size=I(2), alpha=0.75, width=0.2, height=0) +
    ggtitle(paste0(name, " Distribution Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Trait Value") +
    theme_linedraw(base_size=18)
  
  q <- ggplot() +
    geom_bar(data=dists.tbl, aes(x=Population, fill=Trait), position=position_fill(), color="black") +
    ggtitle(paste0(name, " Proportions Across Subtypes"), toupper(COHORT)) +
    xlab("Population") +
    ylab("Proportion") +
    theme_linedraw(base_size=18)
  
  # Add means to continuous variables
  if (!categorical) {
    means <- as.data.frame(dists.tbl %>% group_by(Population) %>% summarise(mean = mean(Trait)))
    rownames(means) <- means$Population
    means <- means[dists.tbl$Population,]
    p <- p + geom_point(data=means, aes(x=Population, y=mean), pch=4, size=I(4))
  }
  
  # Run Single-Factor ANOVA
  dists.tbl$Trait = as.numeric(dists.tbl$Trait)
  res.aov <- aov(Trait ~ Population, data=dists.tbl)
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  
  # If categorical, run a Chi-Square test
  # If not, run a T-Test
  if (categorical) {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      c.tbl <- table(dists.tbl$Trait, dists.tbl$Population)[,c("Total", LETTERS[i])]
      p.vals <- c(p.vals, chisq.test(c.tbl)$p.value)
    }
  } else {
    p.vals <- c()
    for (i in 1:length(levels(pca.cluster.case$Best.partition))) {
      p.vals <- c(p.vals, t.test(subset(dists.tbl, Population=="Total")$Trait, subset(dists.tbl, Population==LETTERS[i])$Trait)$p.value)
    }
  }
  
  p <- p + geom_vline(xintercept=(1+1:length(levels(pca.cluster.case$Best.partition)))[p.vals < 0.05], lty=2, size=1)
  
  # Create bar graph
  if (categorical) {
    q <- q + geom_text(
      aes(x=1+1:length(levels(pca.cluster.case$Best.partition)), y=1, label=ifelse(p.vals < 0.05, "*", "")), 
      vjust=0.2, color="black", position=position_dodge(0.9), size=7.5
    )
    return(list(dist.plot=p, prop.plot=q))
  }
  
  return(p)
}

# Compound submodules by annotation
rosmap.cytokine.signaling <- function(eigengenes) {
  return(svd(eigengenes[,c("DLPFCblue_2", "DLPFCblue_4")])$u[,1])
}
rosmap.leukocyte.activation <- function(eigengenes) {
  return(eigengenes$DLPFCblue_3)
}
rosmap.oligodendrocyte.glial.cells <- function(eigengenes) {
  return(eigengenes$DLPFCbrown_1)
}
rosmap.neuronal <- function(eigengenes) {
  return(svd(eigengenes[,c("DLPFCyellow_1", "DLPFCyellow_2", "DLPFCyellow_3")])$u[,1])
}
rosmap.interleukin.ecm.interactions <- function(eigengenes) {
  return(eigengenes$DLPFCblue_1)
}
rosmap.unfolded.protein <- function(eigengenes) {
  return(svd(eigengenes[,c("DLPFCturquoise_1", "DLPFCturquoise_2")])$u[,1])
}
rosmap.ecm <- function(eigengenes) {
  return(eigengenes$DLPFCbrown_2)
}
rosmap.compound <- function(eigengenes) {
  ret <- data.frame(
    Submodule=colnames(eigengenes)[-(1:3)],
    Annotation=c("Unfolded Protein", "Unfolded Protein", "Interleukin and ECM", "Cytokine Signaling", "Leukocyte Activation", 
                 "Cytokine Signaling", "Oligodendrocyte and Glial Cells", "ECM", "Neuronal", "Neuronal", "Neuronal")
  )
  rownames(ret) <- ret$Submodule
  return(ret)
}
rosmap.sig.genes <- c("APOE", "TREM2", "TYROBP", "BIN2", "APH1B", "AP4M1", 
                      "ABCA7", "MS4A6E", "IL34", "SQSTM1", "SORL1", "GAL3ST4", 
                      "MTMR3", "PTPRB", "APOB", "CLU", "TARDBP", "TMEM106B",
                      "ADAMTS4", "CR1", "BIN1", "INPPD4", "HESX1", "CLNK", 
                      "HS3ST1", "CD2AP", "ZCWPW1",
                      "ECHDC3", "MS4A6A", "PICALM", "SLC24A4", "ADAM10", "APH1B",
                      "KAT8", "SCIMP", "ABI3", "BZRAP1-AS1", "SUZ12P1", "ALPK2",
                      "CD33", "CASS4", "LPL", "CST7", "AXL", "ITGAX", "SPP1", "CD9",
                      "CCL6", "CSF1", "TARDBP", "C9orf72", "GRN", "TBK1", "VCP")

rosmap.cell.astrocytes <- function(eigengenes) svd(eigengenes[,c("DLPFCblue_1", "DLPFCblue_4")])$u[,1]
rosmap.cell.endothelial <- function(eigengenes) eigengenes$DLPFCblue_2
rosmap.cell.microglia <- function(eigengenes) eigengenes$DLPFCblue_3
rosmap.cell.neurons <- function(eigengenes) eigengenes$DLPFCyellow_1
rosmap.cell.oligodendrocytes <- function(eigengenes) eigengenes$DLPFCbrown_1

annot.eigens <- function(eigengenes) {
  mtx = cbind(
    rosmap.cytokine.signaling(eigengenes),
    rosmap.leukocyte.activation(eigengenes),
    rosmap.oligodendrocyte.glial.cells(eigengenes),
    rosmap.neuronal(eigengenes),
    rosmap.interleukin.ecm.interactions(eigengenes),
    rosmap.unfolded.protein(eigengenes),
    rosmap.ecm(eigengenes)
  ) %>% as.data.frame
  rownames(mtx) <- eigengenes$Patient
  colnames(mtx) <- c("Cytokine Signaling", "Leukocyte Activation", "Oligodendrocyte and Glial Cells", "Neuronal", "Interleukin and ECM", "Unfolded Proteins", "ECM")
  return(mtx)
}

cell.eigens <- function(eigengenes) {
  mtx = cbind(
    rosmap.cell.astrocytes(eigengenes),
    rosmap.cell.endothelial(eigengenes),
    rosmap.cell.microglia(eigengenes),
    rosmap.cell.neurons(eigengenes),
    rosmap.cell.oligodendrocytes(eigengenes)
  ) %>% as.data.frame
  rownames(mtx) <- eigengenes$Patient
  colnames(mtx) <- c("Astroctyes", "Endothelial Cells", "Microglial", "Neurons", "Oligodendrocytes")
  return(mtx)
}
