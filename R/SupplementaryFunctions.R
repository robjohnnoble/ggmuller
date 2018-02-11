#' Add branches of length zero to get rid of single nodes in an adjacency matrix
#'
#' Single nodes are those with exactly one daughter. 
#' This function is required by adj_matrix_to_tree, 
#' since valid "phylo" objects cannot contain single nodes. 
#' If pre-existing branches lack lengths then these are set to 1.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#'
#' @return A dataframe comprising the augmented adjacency matrix.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#'
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3), Identity = 2:5)
#' branch_singles(edges1)
#'
#' @export
#' @import dplyr
branch_singles <- function(edges) {
  new_rows <- edges %>% group_by_(~Parent) %>% filter(n() == 1) %>% ungroup()
  num_new_rows <- dim(new_rows)[1]
  if(num_new_rows == 0) return(edges)
  new_rows$Identity <- 1:num_new_rows + max(edges)
  if(is.null(edges$edge.length)) edges$edge.length <- 1
  new_rows$edge.length <- 0
  return(rbind(edges, new_rows))
}

#' Extract an adjacency matrix from a larger data frame
#'
#' @param df Dataframe inclduing column names "Identity", "Parent", and either "Generation" or "Time"
#' @param generation Numeric value of Generation (or Time) at which to determine the adjacency matrix (defaults to final time point)
#'
#' @return A dataframe comprising the adjacency matrix.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_population_df}}
#'
#' @examples
#' edges <- get_edges(example_df)
#' 
#' # extract the adjacency matrix from the data frame:
#' pop_df <- get_population_df(example_df)
#' 
#' # create data frame for plot:
#' Muller_df <- get_Muller_df(edges, pop_df)
#'
#' require(RColorBrewer) # for the palette
#' 
#' # draw plot:
#' num_cols <- length(unique(Muller_df$RelativeFitness)) + 1
#' Muller_df$RelativeFitness <- as.factor(Muller_df$RelativeFitness)
#' Muller_plot(Muller_df, colour_by = "RelativeFitness", 
#'             palette = rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(num_cols)), 
#'             add_legend = TRUE)
#'
#' @export
#' @import dplyr
get_edges <- function(df, generation = NA) {
  if("Time" %in% colnames(df) && !("Generation" %in% colnames(df))) colnames(df)[colnames(df) == "Time"] <- "Generation"
  # check column names:
  if(!("Generation" %in% colnames(df)) | !("Identity" %in% colnames(df)) | !("Parent" %in% colnames(df))) 
    stop("colnames(df) must contain Generation (or Time), Identity and Parent")
  if(is.na(generation)) generation <- max(df$Generation)
  edges <- filter_(df, ~Generation == generation) %>% select_(~Parent, ~Identity)
  edges <- filter_(edges, "Parent != Identity") # remove any row that connects a node to itself
  return(edges)
}

#' Extract population data from a larger data frame
#'
#' @param df Dataframe inclduing column names "Identity", "Parent", and either "Generation" or "Time"
#'
#' @return A dataframe comprising the population dynamics.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_edges}}
#'
#' @examples
#' # extract the adjacency matrix from the data frame:
#' edges <- get_edges(example_df)
#' 
#' # extract the populations (and any other attributes) from the data frame:
#' pop_df <- get_population_df(example_df)
#' 
#' # create data frame for plot:
#' Muller_df <- get_Muller_df(edges, pop_df)
#' 
#' require(RColorBrewer) # for the palette
#' 
#' # draw plot:
#' num_cols <- length(unique(Muller_df$RelativeFitness)) + 1
#' Muller_df$RelativeFitness <- as.factor(Muller_df$RelativeFitness)
#' Muller_plot(Muller_df, colour_by = "RelativeFitness", 
#'             palette = rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(num_cols)), 
#'             add_legend = TRUE)
#'
#' @export
#' @import dplyr
get_population_df <- function(df) {
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(df) && !("Generation" %in% colnames(df))) {
    colnames(df)[colnames(df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # check column names:
  if(!("Generation" %in% colnames(df)) | !("Identity" %in% colnames(df)) | !("Population" %in% colnames(df))) 
    stop("colnames(df) must contain Generation (or Time), Identity and Population")
  
  . <- NULL # avoid check() note
  Population <- NULL # avoid check() note
  
  max_gen <- max(df$Generation)
  max_gen_ids <- filter_(df, ~Generation == max_gen)$Identity
  df <- filter_(df, ~Identity %in% max_gen_ids)
  n <- length(unique(df$Identity))
  master <- data.frame(Generation = rep(unique(df$Generation), each = n),
                       Identity = unique(df$Identity))
  res <- master %>%
    left_join(., df, by = c("Generation", "Identity")) %>%
    mutate(Population = ifelse(Population %in% NA, 0, Population))
  cols <- colnames(df)[!(colnames(df) %in% c("Generation", "Identity", "Population"))]
  for(col in cols) res[, col] <- res[res$Generation == max(res$Generation), col]
  
  # restore original time column name:
  colnames(res)[colnames(res) == "Generation"] <- original_colname
  
  return(res)
}

#' Create a tree object of class "phylo" from an adjacency matrix
#'
#' @param edges Dataframe comprising an adjacency matrix, in which the first column is the parent and the second is the daughter.
#'
#' @return A phylo object.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#'
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' require(ape)
#' tree <- adj_matrix_to_tree(edges1)
#' plot(tree)
#'
#' @export
#' @import dplyr
#' @import ape
adj_matrix_to_tree <- function(edges) {
  # initialise:
  edges <- as.data.frame(edges)
  colnames(edges) <- c("Parent", "Identity")
  if(is.factor(edges$Identity)) {
    edges$Identity <- levels(edges$Identity)[edges$Identity]
    edges$Parent <- levels(edges$Parent)[edges$Parent]
  }
  if(is.character(edges$Identity)) {
    lev_set <- unique(c(edges$Identity, edges$Parent))
    edges$Identity <- factor(edges$Identity, levels = lev_set)
    edges$Parent <- factor(edges$Parent, levels = lev_set)
    edges$Identity <- as.numeric(edges$Identity)
    edges$Parent <- as.numeric(edges$Parent)
  }
  edges <- filter_(edges, "Parent != Identity") # remove any row that connects a node to itself
  edges <- branch_singles(edges) # add branches of length zero to get rid of single nodes
  depth <- 0
  next_rank <- vector()
  n <- 1
  max_depth <- 0
  start_node <- find_start_node(edges)
  path <- start_node
  edges$depth <- NA
  edges$is_tip <- FALSE
  edges[edges$Identity == path[n], "depth"] <- depth
  upped <- FALSE
  num_tips <- 0
  
  # traverse tree, recording depths and ranks:
  repeat {
    if(!upped) repeat { # a downwards move should never follow an upwards move
      n <- n + 1
      path[n] <- move_down(edges, path[n - 1])
      if(path[n] != path[n - 1]) {
        depth <- depth + 1
        if(depth > max_depth) {
          next_rank[depth] <- 0
          max_depth <- depth
        }
        next_rank[depth] <- next_rank[depth] + 1
      }
      else {
        edges[edges$Identity == path[n], "is_tip"] <- TRUE
        num_tips <- num_tips + 1
      }
      edges[edges$Identity == path[n], "depth"] <- depth
      edges[edges$Identity == path[n], "rank_at_depth"] <- next_rank[depth]
      upped <- FALSE
      if(path[n] == path[n - 1]) break
    }
    if(move_right(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_right(edges, path[n - 1])
      if(path[n] != path[n - 1]) next_rank[depth] <- next_rank[depth] + 1
      edges[edges$Identity == path[n], "rank_at_depth"] <- next_rank[depth]
      edges[edges$Identity == path[n], "depth"] <- depth
      upped <- FALSE
    } else if(move_up(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_up(edges, path[n - 1])
      if(path[n] != path[n - 1]) depth <- depth - 1
      edges[edges$Identity == path[n], "depth"] <- depth
      upped <- TRUE
    }
    if(path[n] == path[1]) break
    if(n > 1E6) stop("Error: stuck in a loop")
    if(max(table(path) > 2)) stop("Error: adjacency matrix seems to include loops.")
  }
  if(length(path) != 2 * dim(edges)[1] + 2) stop("Error: adjacency matrix seems to be bipartite.")
  
  # relabel nodes to conform with "phylo" standard:
  edges <- rbind(filter_(edges, ~is_tip) %>% arrange_(~-depth, ~rank_at_depth), filter_(edges, ~is_tip == FALSE) %>% arrange_(~depth, ~rank_at_depth))
  if(length(unique(edges$Identity)) > num_tips) edges$New_Identity <- c(1:num_tips, (num_tips + 2):(length(unique(edges$Identity)) + 1))
  else edges$New_Identity <- c(1:num_tips)
  edges$New_Parent <- NA
  for(i in 1:dim(edges)[1]) { # to do: replace this for loop with something more efficient!
    if(edges[i, "Parent"] == start_node) edges[i, "New_Parent"] <- num_tips + 1
    else edges[i, "New_Parent"] <- edges[edges$Identity == edges[i, "Parent"], "New_Identity"]
  }
  
  # create phylo object:
  tree <- list()
  tree$edge.length <- edges$edge.length
  edges <- select_(edges, ~New_Parent, ~New_Identity)
  colnames(edges) <- NULL
  rownames(edges) <- NULL
  tree$edge <- as.matrix(edges)
  tree$edge <- cbind(as.integer(tree$edge[,1]), as.integer(tree$edge[,2]))
  tree$Nnode <- as.integer(max(edges) - num_tips)
  tree$tip.label <- as.character(rep(NA, num_tips))
  class(tree) <- "phylo"
  reorder.phylo(tree, order = "cladewise")
  attr(tree, "order") <- "cladewise"
  
  return(tree)
}
