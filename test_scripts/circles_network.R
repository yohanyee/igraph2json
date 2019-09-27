library(tidyverse)
library(igraph)
library(packcircles)
library(ggrepel)
library(glue)

file.remove(list.files(path = "/tmp", pattern = "fbp_", full.names = T))

# Set params
edge_buffer_fraction <- 0.1
center_exponent <- 1
radial_exponent <- 6
mult_scale <- 0.1
internal_attractive_scale <- 5
internal_attractive_exponent <- 4
internal_repulsive_scale <- 0.2
internal_repulsive_exponent <- 1
external_attractive_scale <- 5
external_attractive_exponent <- 4
force_step <- 0.0005

# Create initial data
modules <- cluster_walktrap(g)

num_modules <- length(unique(modules$membership))

nodes <- data.frame(name=modules$names, module=factor(modules$membership), index=1:length(modules$membership)) 
edges <- as_long_data_frame(g) %>%
  select(from, to, weight)

# Initial module ordering 
module_adjacency <- matrix(nrow=num_modules, ncol=num_modules)
for (i in 1:num_modules) {
  for (j in 1:num_modules) {
    module_adjacency[i,j] <- mean(as.matrix(as_adj(g, attr = "weight")[modules$membership==i, modules$membership==j]))
  }
}
module_adjacency <- (module_adjacency + t(module_adjacency))/2
module_order <- hclust(dist(module_adjacency), method="complete")$order
module_sizes <- table(modules$membership)

# Initial module packing
module_packing <- circleProgressiveLayout(module_sizes[module_order]) %>%
  mutate(module=factor(module_order)) %>%
  rename(module_center_x=x,
         module_center_y=y,
         module_radius=radius)

# Functions

get_central_force <- function(x, y, center_x, center_y, radius, center_exponent=1, radial_exponent=1, mult_scale=0.01) {
  
  a <- center_exponent
  b <- radial_exponent
  m <- mult_scale
  x_diff <- x - center_x
  y_diff <- y - center_y
  
  d <- sqrt((x_diff)^2 + (y_diff)^2)
  scaled_d <- d/radius
  if (scaled_d < 1) {
    force_magnitude <- (m*(a*(scaled_d-1) + b*scaled_d)/((scaled_d-1)^(b+1)*scaled_d^(a+1)))
  } else {
    force_magnitude <- -100
  }
  
  force_x <- force_magnitude * x_diff/sqrt(x_diff^2 + y_diff^2)
  force_y <- force_magnitude * y_diff/sqrt(x_diff^2 + y_diff^2)
  
  return(c(force_x, force_y))
}


get_node_internal_force <- function(x, y, x_internal, y_internal, radius, weight, attractive_scale=5, attractive_exponent=4, repulsive_scale=0.01, repulsive_exponent=1) {
  
  aa <- attractive_scale
  ba <- attractive_exponent
  ar <- repulsive_scale
  br <- repulsive_exponent
  
  x_diff <- x - x_internal
  y_diff <- y - y_internal
  
  d <- sqrt((x-x_internal)^2 + (y-y_internal)^2)
  scaled_d <- d/radius
  
  force_magnitude_repulsive <- ar*br/scaled_d^(br+1)
  force_magnitude_attractive <- -aa*ba*weight^(ba-1)
  
  force_repulsive_x <- force_magnitude_repulsive * x_diff/sqrt(x_diff^2 + y_diff^2)
  force_repulsive_y <- force_magnitude_repulsive * y_diff/sqrt(x_diff^2 + y_diff^2)
  force_attractive_x <- force_magnitude_attractive * x_diff/sqrt(x_diff^2 + y_diff^2)
  force_attractive_y <- force_magnitude_attractive * y_diff/sqrt(x_diff^2 + y_diff^2)
  
  return(c(force_attractive_x + force_repulsive_x, force_attractive_y + force_repulsive_y))
}

get_node_external_force <- function(x, y, x_external, y_external, weight, attractive_scale=5, attractive_exponent=4) {
  aa <- attractive_scale
  ba <- attractive_exponent

  x_diff <- x - x_external
  y_diff <- y - y_external
  
  force_magnitude_attractive <- -aa*ba*weight^(ba-1)

  force_attractive_x <- force_magnitude_attractive * x_diff/sqrt(x_diff^2 + y_diff^2)
  force_attractive_y <- force_magnitude_attractive * y_diff/sqrt(x_diff^2 + y_diff^2)
  
  return(c(force_attractive_x, force_attractive_y))
}


# State
state <- nodes %>% 
  left_join(module_packing, by="module") %>%
  mutate(x=rnorm(length(module_center_x), mean=module_center_x, sd=module_radius/10),
         y=rnorm(length(module_center_y), mean=module_center_y, sd=module_radius/10))




##########
for (iter in 1:250) {
  
  total_movement <- 0
  
  # Loop over each node
  for (i in 1:nrow(state)) {
    
    previous_state <- state
    
    edges_iter <- edges %>%
      select(from, to, weight) %>%
      left_join(state %>% 
                  select(name, module, index, x, y) %>%
                  rename(from_name=name,
                         from_module=module,
                         from=index,
                         from_x=x,
                         from_y=y), 
                by="from") %>%
      left_join(state %>% 
                  select(name, module, index, x, y) %>%
                  rename(to_name=name,
                         to_module=module,
                         to=index,
                         to_x=x,
                         to_y=y), 
                by="to")
    
    
    # Central force relative to module circle
    central_forces <- get_central_force(x = state$x[i], 
                      y = state$y[i], 
                      center_x = state$module_center_x[i],
                      center_y = state$module_center_y[i], 
                      radius = state$module_radius[i], 
                      center_exponent = center_exponent, radial_exponent = radial_exponent, mult_scale = mult_scale) %>% 
      `names<-`(., c("fx", "fy"))
    
    # Node force loop over connections
    edges_state <- edges_iter %>%
      filter(from==i, from_module==to_module) 
    
    node_internal_forces_lowerhalf <- 1:nrow(edges_state) %>%
      map_dfr(function(i) {
        f <- get_node_internal_force(x = edges_state$from_x[i], 
                                y = edges_state$from_y[i], 
                                x_internal = edges_state$to_x[i], 
                                y_internal = edges_state$to_y[i], 
                                weight = edges_state$weight[i], 
                                radius = state$module_radius, 
                                attractive_scale = internal_attractive_scale, 
                                attractive_exponent = internal_attractive_exponent, 
                                repulsive_scale = internal_repulsive_scale, 
                                repulsive_exponent = internal_repulsive_exponent)
        return(data.frame(fx=f[1], fy=f[2]))
      }) %>%
      colMeans()
    
    edges_state <- edges_iter %>%
      filter(to==i, from_module==to_module) 
    
    node_internal_forces_upperhalf <- 1:nrow(edges_state) %>%
      map_dfr(function(i) {
        f <- get_node_internal_force(x = edges_state$to_x[i], 
                                     y = edges_state$to_y[i], 
                                     x_internal = edges_state$from_x[i], 
                                     y_internal = edges_state$from_y[i], 
                                     weight = edges_state$weight[i], 
                                     radius = state$module_radius, 
                                     attractive_scale = internal_attractive_scale, 
                                     attractive_exponent = internal_attractive_exponent, 
                                     repulsive_scale = internal_repulsive_scale, 
                                     repulsive_exponent = internal_repulsive_exponent)
        return(data.frame(fx=f[1], fy=f[2]))
      }) %>%
      colMeans()
    
    
    # Node force loop over connections
    edges_state <- edges_iter %>%
      filter(from==i, from_module != to_module) 
    
    node_external_forces_lowerhalf <- 1:nrow(edges_state) %>%
      map_dfr(function(i) {
        f <- get_node_external_force(x = edges_state$from_x[i], 
                                     y = edges_state$from_y[i], 
                                     x_external = edges_state$to_x[i], 
                                     y_external = edges_state$to_y[i], 
                                     weight = edges_state$weight[i], 
                                     attractive_scale = external_attractive_scale, 
                                     attractive_exponent = external_attractive_exponent)
        return(data.frame(fx=f[1], fy=f[2]))
      }) %>%
      colMeans()
    
    edges_state <- edges_iter %>%
      filter(to==i, from_module != to_module) 
    
    node_external_forces_upperhalf <- 1:nrow(edges_state) %>%
      map_dfr(function(i) {
        f <- get_node_external_force(x = edges_state$to_x[i], 
                                     y = edges_state$to_y[i], 
                                     x_external = edges_state$from_x[i], 
                                     y_external = edges_state$from_y[i], 
                                     weight = edges_state$weight[i], 
                                     attractive_scale = external_attractive_scale, 
                                     attractive_exponent = external_attractive_exponent)
        return(data.frame(fx=f[1], fy=f[2]))
      }) %>%
      colMeans()
    
    
    total_forces <- rbind(central_forces,
          node_internal_forces_lowerhalf, 
          node_internal_forces_upperhalf,
          node_external_forces_lowerhalf,
          node_external_forces_upperhalf) %>%
      colSums(na.rm = T)
    
    this_step <- force_step*c(total_forces[1], total_forces[2]) 
    total_movement <- total_movement + this_step[1]^2 + this_step[2]^2
    
    state$x[i] <- state$x[i] + sign(this_step[1])*min(abs(this_step[1]), state$module_radius[i]/10)
    state$y[i] <- state$y[i] + sign(this_step[1])*min(abs(this_step[2]), state$module_radius[i]/10)
    
    if (any(is.nan(state$x)) | any(is.nan(state$y))) {
      print(i)
    }
  }
  
  mean_abs_x_movement <- mean(abs(state$x - previous_state$x))
  mean_abs_y_movement <- mean(abs(state$y - previous_state$y))
  max_point_movement <- max(abs(state$x - previous_state$x), abs(state$y - previous_state$y))
  
  cat(glue("Completed iteration {iter}. Max movement is {max_point_movement}\n", .trim=F))
  cat(glue("Mean abs(delta x) = {mean_abs_x_movement}\n", .trim=F))
  cat(glue("Mean abs(delta y) = {mean_abs_y_movement}\n", .trim=F))
  cat("\n")
  
  plt_module_polygon <- circleLayoutVertices(module_packing, npoints=50, idcol = "module", xysizecols = c(1,2,3)) %>%
    rename(module=id) %>%
    mutate(module=factor(module))
  
  plt_edges <- edges %>% 
    select(from, to, weight) %>%
    left_join(state %>% 
                select(name, module, index, x, y) %>%
                rename(from_name=name,
                       from_module=module,
                       from=index,
                       from_x=x,
                       from_y=y), 
              by="from") %>%
    left_join(state %>% 
                select(name, module, index, x, y) %>%
                rename(to_name=name,
                       to_module=module,
                       to=index,
                       to_x=x,
                       to_y=y), 
              by="to")
  
  png(glue("/tmp/fbp_{sprintf('%03d', iter)}.png"))
  print(  # Plot example
    # Plot

    ggplot() + 
      geom_polygon(data = plt_module_polygon,
                   aes(x=x, y=y, group = module, colour=module), 
                   fill = "white", 
                   alpha = 0.7, 
                   show.legend = FALSE,
                   linetype='dashed') +
      geom_segment(data = plt_edges %>%
                     filter(from_module==to_module), 
                   aes(x=from_x, y=from_y, xend=to_x, yend=to_y, color=from_module, alpha=weight),
                   size=0.5,
                   lineend = "round",
                   linejoin = "mitre",
                   arrow = arrow(length = unit(0.01, "npc"))) + 
      geom_curve(data = plt_edges %>%
                   filter(from_module != to_module), 
                 aes(x=from_x, y=from_y, xend=to_x, yend=to_y, color=from_module, alpha=weight),
                 size=0.5,
                 lineend = "round",
                 arrow = arrow(length = unit(0.01, "npc"))) + 
      geom_point(data = state,
                 aes(x=x, y=y, fill=module),
                 color="black", 
                 shape=21,
                 size=2,
                 stroke=1) +
      coord_fixed(ratio=1) +
      guides(color=F, alpha=F, fill=F) +
      theme_void()
    )
  dev.off()
}





##########
# Plot
plt_module_polygon <- circleLayoutVertices(module_packing, npoints=50, idcol = "module", xysizecols = c(1,2,3)) %>%
  rename(module=id) %>%
  mutate(module=factor(module))

plt_edges <- edges %>% 
  select(from, to, weight) %>%
  left_join(state %>% 
              select(name, module, index, x, y) %>%
              rename(from_name=name,
                     from_module=module,
                     from=index,
                     from_x=x,
                     from_y=y), 
            by="from") %>%
  left_join(state %>% 
              select(name, module, index, x, y) %>%
              rename(to_name=name,
                     to_module=module,
                     to=index,
                     to_x=x,
                     to_y=y), 
            by="to")

ggplot() + 
  geom_polygon(data = plt_module_polygon,
               aes(x=x, y=y, group = module, colour=module), 
               fill = "white", 
               alpha = 0.7, 
               show.legend = FALSE,
               linetype='dashed') +
  geom_segment(data = plt_edges %>%
                 filter(from_module==to_module), 
               aes(x=from_x, y=from_y, xend=to_x, yend=to_y, color=from_module, alpha=weight),
               size=0.5,
               lineend = "round",
               linejoin = "mitre",
               arrow = arrow(length = unit(0.01, "npc"))) + 
  geom_curve(data = plt_edges %>%
               filter(from_module != to_module), 
             aes(x=from_x, y=from_y, xend=to_x, yend=to_y, color=from_module, alpha=weight),
             size=0.5,
             lineend = "round",
             arrow = arrow(length = unit(0.01, "npc"))) + 
  geom_label_repel(data = state %>%
                     select(name, module, index, x, y), 
                   aes(label=name, x=x, y=y),
                   segment.alpha = 0.8,
                   force = 10) +
  geom_point(data = state,
             aes(x=x, y=y, fill=module),
             color="black", 
             shape=21,
             size=2,
             stroke=1) +
  coord_fixed(ratio=1) +
  guides(color=F, alpha=F, fill=F) +
  theme_void()


# To implement
# labeling of modules
# -manual
# -via enrichment
# -curved vs straight connectors (separate for inner vs outer connections)
# Positioning of modules
# -manual


