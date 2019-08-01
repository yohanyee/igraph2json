library(tidyverse)
library(magrittr)
library(igraph)
library(corrplot)

#######################
# Get example data ----
#######################

gene_symbols <- c("Fmr1")
min_score <- 700 # Out of 1000

# Get STRINGids for an example gene
api_call <- paste("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(gene_symbols, collapse = "%0D"), "&species=10090&limit=1", sep="")
string_ids <- read_tsv(api_call) %>% select(queryIndex, stringId, preferredName)

# Remove organism ID
genes <- string_ids$stringId %>% strsplit("\\.") %>% map_chr(`[`, 2)

# Get interaction partners (nodes)
api_call <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=", paste(genes, collapse = "%0D"), "&species=10090&required_score=", min_score, sep="")
interaction_data <- read_tsv(api_call)

# Get network between interaction partners
genes <- unique(c(interaction_data$stringId_A, interaction_data$stringId_B))
api_call <- paste("https://string-db.org/api/tsv/network?identifiers=", paste(genes, collapse = "%0D"), "&species=10090&required_score=", min_score, sep="")
network_data <- read_tsv(api_call)

# Construct graph
node_df <- data.frame(id=unique(c(network_data$preferredName_A, network_data$preferredName_B)))
edge_df <- network_data %>% 
  select(preferredName_A, preferredName_B, score) %>%
  rename(source=preferredName_A, target=preferredName_B)

g <- graph_from_data_frame(edge_df, vertices = node_df)
E(g)$weight <- get.edge.attribute(g)$score
attributes(g)

# Detect communities
communities <- cluster_walktrap(g)
node_df <- node_df %>% inner_join(data.frame(id=names(membership(communities)), 
                                  community=as.integer(membership(communities))),
                       by="id")

# Also create adjacency matrix
adj <- as_adj(g,attr = "score")
corrplot(as.matrix(adj), is.corr=F, method="color")

# Display
plot(g, layout=layout.fruchterman.reingold(g), edge.width=(E(g)/100)^2, edge.arrow.mode=0, vertex.label.family="sansserif")
plot(communities, g, layout=layout.fruchterman.reingold(g), edge.width=(E(g)/100)^2, edge.arrow.mode=0, vertex.label.family="sansserif")

# Also save as png
png("graph_example_STRINGdb_igraph.png", width = 1024, height = 786)
print(plot(g, layout=layout.fruchterman.reingold(g), edge.width=(E(g)/100)^2, edge.arrow.mode=0, vertex.label.family="sansserif"))
dev.off()

png("graph_example_STRINGdb_igraph_with_communities.png", width = 1024, height = 786)
print(plot(communities, g, layout=layout.fruchterman.reingold(g), edge.width=(E(g)/100)^2, edge.arrow.mode=0, vertex.label.family="sansserif"))
dev.off()

# Write out
source("R/convert.R")

# If creating plots directly from dataframes that contain node and edge information
create_plot_from_df(node_df = node_df %>% mutate(name=id), edge_df = edge_df, templatedir="template", outdir="plot_from_df", node_size_column = NA, node_size = 5)

# If creating plots from igraph object
create_plot_from_igraph(g, templatedir="template", outdir="plot_from_igraph", node_size_column = NA, node_size = 5)

# If creating plots from adjacency matrix
create_plot_from_adj(adj, templatedir="template", outdir="plot_from_adj", node_size_column = NA, node_size = 5)
