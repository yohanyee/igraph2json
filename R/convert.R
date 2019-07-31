library(tidyverse)
library(magrittr)
library(igraph)
library(rjson)

edges <- E(g)
nodes <- V(g)

# Adjacency matrix to graph

# Graph to dataframe

# Dataframe to json
# node_df: must have columns "id" and optionally "community"
# edge_df: must have columns "source" and "target"
df_to_json <- function(node_df, edge_df) {
  # Initialize list object
  jsobj <- list()

  # Add nodes
  jsobj$nodes <- list()
  for (i in 1:nrow(node_df)) {
    jsobj$nodes[[i]] <- as.list(node_df[i,])
    
  }
  
  # Add edges
  jsobj$links <- list()
  for (i in 1:nrow(edge_df)) {
    jsobj$links[[i]] <- as.list(edge_df[i,])
  }
  
  return(toJSON(jsobj))
}

#####################
# Write out data ----
#####################

# Write json
write_json <- function(x, path) {
  write_file(x, path)  
}

# Write tempate
write_template <- function(template_path, output_path, jsonfile, node_size_column=NULL, link_width_column=NULL, node_size=5, link_width=1) {
  template_string <- read_file(template_path)
  
  template_string <- gsub("%%FILE%%", jsonfile, template_string)
  
  if (!is.null(node_size_column)) {
    template_string <- gsub("%%GRAPH_NODE_SORT%%", paste0("graph.nodes.sort(function (a, b) { return b.", node_size_column, " - a.", node_size_column, "; });"), template_string)
    template_string <- gsub("%%NODE_RADIUS_DOMAIN%%", paste0("nodeRadius.domain([graph.nodes[graph.nodes.length-1].", node_size_column, ", graph.nodes[0].", node_size_column, "]);"), template_string)
    template_string <- gsub("%%PLOT_NODE_ATTR_NODE_RADIUS%%", paste0(".attr('r', function (d) { return nodeRadius(d.", node_size_column, "); })"), template_string)
  } else {
    template_string <- gsub("%%GRAPH_NODE_SORT%%", "", template_string)
    template_string <- gsub("%%NODE_RADIUS_DOMAIN%%", "", template_string)
    template_string <- gsub("%%PLOT_NODE_ATTR_NODE_RADIUS%%", paste0(".attr('r', function (d) { return ", node_size, "; })"), template_string)
  }
  
  if (!is.null(link_width_column)) {
    template_string <- gsub("%%LINK_WIDTH_DOMAIN%%", paste0("linkWidth.domain(d3.extent(graph.links, function (d) { return d.", link_width_column, "; }));"), template_string)
    template_string <- gsub("%%PLOT_LINK_ATTR_STROKE_WIDTH%%", paste0(".attr('stroke-width', function (d) { return linkWidth(d.", link_width_column, "); });"), template_string)
  } else {
    template_string <- gsub("%%LINK_WIDTH_DOMAIN%%", "", template_string)
    template_string <- gsub("%%PLOT_LINK_ATTR_STROKE_WIDTH%%", paste0(".attr('stroke-width', function (d) { return ", link_width, "; });"), template_string)
    
  }
  
  write_file(template_string, output_path)
}

# Write dataframe
write_df <- function(templatedir, outdir, node_df, edge_df, node_size_column=NULL, link_width_column=NULL, node_size=5, link_width=1) {
  
  jsobj <- df_to_json(node_df, edge_df)
  
  if (is.null(node_size_column) & ncol(node_df) >= 2) {
    cnames <- colnames(node_df)
    if ("name" %in% cnames) {
      cnames <- cnames[-c(which(cnames=="name"))]
    }
    node_size_column <- cnames[2]
  } else {
    if (is.na(node_size_column)) {
      node_size_column <- NULL
    }
  }
  
  if (is.null(link_width_column) & ncol(edge_df) >= 3) {
    link_width_column <- colnames(edge_df)[3]
  } else {
    if (is.na(link_width_column)) {
      link_width_column <- NULL
    }
  }
  
  if (templatedir != outdir) {
    dir.create(outdir, showWarnings = F, recursive = T)  
    
    write_template(file.path(templatedir, "template.html"), 
                   file.path(outdir, "index.html"), 
                   "graph.json", 
                   node_size_column, 
                   link_width_column,
                   node_size,
                   link_width)
    write_json(jsobj, 
               file.path(outdir, "graph.json"))
    
    file.copy(file.path(templatedir, "forceInABox.js"), file.path(outdir, "forceInABox.js"), overwrite = T)  
  }
  
}



