# Goal
# * Develop functions to plot bayesian networks produced by bnlearn

# Load libraries
library(bnlearn)
library(igraph)
library(dplyr)

# Takes a vector of intensity value (between 0 and 1) 
# and converts it to a hex representation between white and black
intensity.to.color <- function(intensity) {
  
  # Interpolate between blue and red and convert to hex values
  interpolated <- colorRamp(c("#3300FF", "#FF0033"))(intensity) / 255
  hexes <- character(nrow(interpolated))
  # Store interpolated values as hexes
  for (i in 1:nrow(interpolated)) {
    hexes[i] <- rgb(
      interpolated[i,1],
      interpolated[i,2],
      interpolated[i,3]
    )
  }
  return(hexes)
}

# Takes an object produced by boot.strength and
# produces a vector of unique nodes for iGraph
bn.to.labels <- function(bn) {
  
  return(unique(bn$from))
}

# Takes an object produced by boot.strength and 
# produces an object that can be used for iGraph
bn.to.edges <- function(bn, sig.labels) {
  
  label.set <- bn.to.labels(bn) # Retrieve set of labels
  edge.set <- c(rbind(bn$from, bn$to)) %>% as.factor() # Retrieve unique set of edges
  levels(edge.set) <- label.set # Make sure labels and edges are consistent
  edge.set <- as.numeric(edge.set) # Coerce into a numeric form
  
  edges <- list() # List of edges selected
  edges.width <- list() # Widths of edges
  edges.arrow.width <- list() # Widths of arrows
  edges.count <- 0 # Iterator
  rows.chosen <- c() # For blue-red coloring of edges
  
  for (i in 1:nrow(bn)) {
    
    if (bn[i,]$direction >= 0.5 & bn[i,]$strength >= 0.95) {
      
      edges.count <- edges.count + 1
      
      rows.chosen <- c(rows.chosen, i)
      edges[[edges.count]] <- edge.set[(2*i-1):(2*i)] # Add edge into set
      
      # Check if the arrow is between two significant nodes
      # If it is, increase weight of arrow to show significant subnetwork
      if (bn[i,]$from %in% sig.labels & bn[i,]$to %in% sig.labels) {
        edges.width[[edges.count]] <- 6
        edges.arrow.width[[edges.count]] <- 2
      } else {
        edges.width[[edges.count]] <- 2
        edges.arrow.width[[edges.count]] <- 1
      }
    }
  }
  edges <- unlist(edges)
  edges.width <- unlist(edges.width)
  edges.arrow.width <- unlist(edges.arrow.width)
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))} # Define function to normalize between range
  
  # Create object and store edge properties
  edge.obj <- list(
    edges=edges,
    # edges.color=intensity.to.color(range01(bn[rows.chosen,]$strength)), # For blue-red coloring of edges
    edges.color="#FF0033",
    edges.width=edges.width,
    edges.arrow.width=edges.arrow.width
  )
  
  return(edge.obj)
}

# Takes an object returned by boost.strength and
# determines the colors
bn.to.node.colors <- function(bn) {
  
  label.set <- bn.to.labels(bn)
  color.set <- rep("#FFFFFF", length(label.set))
  color.set[grep("turquoise", label.set)] <- "#40E0D0"
  color.set[grep("blue", label.set)] <- "#ADD8E6"
  color.set[grep("yellow", label.set)] <- "#F8DE7E"
  color.set[grep("green", label.set)] <- "#98FB98"
  color.set[grep("brown", label.set)] <- "#E08B3E"
  
  return(color.set)
}

# Takes the object returned by boot.strength and creates
# a graph representing the bayesian network
plot.bn <- function(bn, sig.labels=c(), v.color="black", v.cex=2, v.font=2, v.size=30, e.arrow.size=4, v.colors=c()) {
  
  # Retrieve lables and edges
  bn.labels <- bn.to.labels(bn)
  bn.edges <- bn.to.edges(bn, sig.labels)
  
  # Create empty graph with labeled nodes
  g <- graph.empty(n=length(bn.labels), directed=T)
  
  # Add all edges 
  g <- add.edges(g, bn.edges$edges)
  
  # Set vertex properties
  V(g)$label <- bn.labels
  V(g)$color <- bn.to.node.colors(bn)
  V(g)$label.color <- v.color
  V(g)$label.cex <- v.cex
  V(g)$label.font <- v.font
  V(g)$size <- v.size
  
  if (length(v.colors) != 0) {
    V(g)$color <- v.colors
  }
  
  # Set edge properties
  E(g)$color <- bn.edges$edges.color
  E(g)$arrow.size <- e.arrow.size
  E(g)$width <- bn.edges$edges.width
  
  return(plot(g, layout=layout.auto))
}