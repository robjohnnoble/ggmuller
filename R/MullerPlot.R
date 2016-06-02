#' Move to daughter in adjacency matrix
#'
#' Returns the smallest Identity value among the set of daughters. 
#' When parent has no daughters, returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param identity Integer specifying parent whose daughter is to be found
#'
#' @return An integer specifying the daughter's Identity.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{move_up}} \code{\link{move_right}}
#'
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' move_down(edges1, 3)
#'
#' @export
#' @import dplyr
move_down <- function(edges, identity) {
  daughters <- filter(edges, Parent == identity)$Identity
  if(length(daughters) == 0) return(identity) # if it is not a parent then don't move
  return(min(daughters))
}

#' Move to sibling in adjacency matrix
#'
#' Returns the next smallest Identity value among the set of siblings. 
#' When there is no such sibling, returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param identity Integer specifying parent whose daughter is to be found
#'
#' @return An integer specifying the sibling's Identity.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{move_up}} \code{\link{move_down}}
#'
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' move_right(edges1, 3)
#'
#' @export
#' @import dplyr
move_right <- function(edges, identity) {
  parent <- filter(edges, Identity == identity)$Parent
  if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
  siblings <- filter(edges, Parent == parent & Identity > identity)$Identity
  if(length(siblings) == 0) return(identity) # if it is the initial genotype then don't move
  return(min(siblings))
}

#' Move to parent in adjacency matrix
#'
#' Returns the corresponding Parent value. 
#' When there is no parent (i.e. at the top of the tree), returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param identity Integer specifying daughter whose parent is to be found
#'
#' @return An integer specifying the Parent.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{move_down}} \code{\link{move_right}}
#' 
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' move_up(edges1, 3)
#'
#' @export
#' @import dplyr
move_up <- function(edges, identity) {
  parent <- filter(edges, Identity == identity)$Parent
  if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
  return(as.numeric(parent))
}

#' Move to top of adjacency matrix
#'
#' Returns the Parent value of the common ancestor.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#'
#' @return An integer specifying the Parent that is the common ancestor.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#'
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' find_start_node(edges1)
#'
#' @export
#' @import dplyr
find_start_node <- function(edges) {
  start <- min(edges$Parent) # reasonable guess
  repeat {
    if(move_up(edges, start) == start) break
    start <- move_up(edges, start)
  }
  return(start)
}

#' Record a path through all nodes of an adjacency matrix
#'
#' Nodes are traversed in the order that they should be stacked in a Muller plot. 
#' Each node appears exactly twice.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#'
#' @return A vector specifying the path.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' 
#' @examples
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
#' get_path(edges1)
#'
#' @export
#' @import dplyr
get_path <- function(edges) {
  n <- 1
  path <- find_start_node(edges)
  upped <- FALSE
  repeat {
    if(!upped) repeat { # a downwards move should never follow an upwards move
      n <- n + 1
      path[n] <- move_down(edges, path[n - 1])
      upped <- FALSE
      if(path[n] == path[n - 1]) break
    }
    if(move_right(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_right(edges, path[n - 1])
      upped <- FALSE
    } else if(move_up(edges, path[n]) != path[n]) {
      n <- n + 1
      path[n] <- move_up(edges, path[n - 1])
      upped <- TRUE
    }
    if(path[n] == path[1]) break
    if(n > 1E6) return("Error: stuck in a loop") # to do: add more sophisticated error checking for loops, bipartite graphs, etc.
    if(max(table(path) > 2)) return("Error: adjacency matrix seems to include loops.")
  }
  if(length(path) != 2 * dim(edges)[1] + 2) return("Error: adjacency matrix seems to be bipartite.")
  return(path)
}

#' Reorder a Muller plot dataframe by a vector
#'
#' @param df Dataframe with column names "Generation", "Parent" and "Identity", in which each Identity appears exactly twice
#' @param vector Vector of Identity values
#'
#' @return The reordered dataframe.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_path}}
#'
#' @examples
#' df <- data.frame(Generation = c(rep(0, 6), rep(1, 6)), 
#'  Identity = rep(1:6,2), Population = c(1, rep(0, 5), 10, rep(1, 5)))
#' df <- rbind(df, df) # duplicate rows
#' require(dplyr)
#' df <- arrange(df, Generation) # put in chronological order
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6) # adjacency matrix
#' path <- get_path(edges1) # path through the adjacency matrix
#' reorder_by_vector(df, path)
#'
#' @export
#' @import dplyr
reorder_by_vector <- function(df, vector) {
  vector <- group_by(as.data.frame(vector), vector) %>% 
    mutate(Unique_id = c(vector[1], paste0(vector[1], "a")))
  gens <- unique(df$Generation)
  n_gens <- length(gens)
  n_ids <- dim(df)[1] / n_gens
  df <- group_by(df, Generation, Identity) %>% 
    mutate(Unique_id = c(paste0(Identity[1], "_", Generation[1]), paste0(Identity[1], "a", "_", Generation[1]))) %>% 
    ungroup
  vector <- rep(vector$Unique_id, n_gens)
  vector <- paste0(vector, "_", rep(gens, each = n_ids))
  df <- df[match(vector, df$Unique_id), ]
  df <- group_by(df, Generation, Identity) %>% 
    mutate(Group_id = c(Identity[1], paste0(Identity[1], "a"))) %>% 
    ungroup
  return(df)
}

#' Create a data frame from which to create a Muller plot
#'
#' @param edges Dataframe comprising an adjacency matrix, or tree of class "phylo"
#' @param pop_df Dataframe with column names "Generation", "Identity" and "Population"
#' @param threshold Numeric threshold; genotypes that never become more abundant than this threshold are omitted
#' @param add_zeroes Logical whether to include rows with Population = 0
#'
#' @return A dataframe that can be used as input in Muller_plot.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{Muller_plot}}
#'
#' @examples
#' # by default, all genotypes are included, 
#' # but one can choose to omit genotypes with max frequency < threshold:
#' Muller_df <- get_Muller_df(example_edges, example_pop_df, threshold = 0.005)
#'
#' # one can also choose to include rows with Population = 0:
#' Muller_df <- get_Muller_df(example_edges, example_pop_df, add_zeroes = TRUE, threshold = 0.005)
#'
#' @export
#' @import dplyr
get_Muller_df <- function(edges, pop_df, add_zeroes = FALSE, threshold = 0) {
  # check/set column names:
  if(!("Generation" %in% colnames(pop_df)) | !("Identity" %in% colnames(pop_df)) | !("Generation" %in% colnames(pop_df))) 
    return("colnames(pop_df) must contain Generation, Identity and Population")
  if(class(edges) == "phylo") edges <- edges$edge
  colnames(edges) <- c("Parent", "Identity")
    
  # add semi-frequencies:
  pop_df <- pop_df %>% group_by(Generation) %>% 
    mutate(Frequency = (Population / sum(Population)) / 2) %>% 
    ungroup
  
  # duplicate rows:
  Muller_df <- rbind(pop_df, pop_df)
  Muller_df <- arrange(Muller_df, Generation)
  
  # get the path:
  path <- get_path(edges)
  path <- rev(path) # apparently, the convention for Muller plots to have earliest-arriving genotypes plotted nearest the top
  
  # rearrange the population data according to the path:
  Muller_df <- reorder_by_vector(Muller_df, path)
  
  # optionally remove rows with Population = 0:
  if(!add_zeroes) Muller_df <- filter(Muller_df, Population > 0)
  
  # optionally remove rare genotypes, and recalculate frequencies:
  if(threshold > 0) {
    Muller_df <- Muller_df %>% group_by(Identity) %>% 
      filter(max(Frequency) >= threshold)
    Muller_df <- Muller_df %>% group_by(Generation) %>% 
      mutate(Frequency = Population / sum(Population)) %>% 
      ungroup
  }
  
  return(Muller_df)
}

#' Draw a Muller plot using ggplot2
#'
#' @param Muller_df Dataframe created by get_Muller_df
#' @param colour_by Character containing name of column by which to colour the plot
#' @param palette List of colours
#' @param add_legend Logical whether to show legend
#'
#' @return None
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_Muller_df}}
#'
#' @examples
#' Muller_df <- get_Muller_df(example_edges, example_pop_df, threshold = 0.005)
#' Muller_plot(Muller_df)
#'
#' @export
#' @import dplyr
#' @import ggplot2
Muller_plot <- function(Muller_df, colour_by = NA, palette = NA, add_legend = FALSE) {
  if(is.na(palette[1])) {
    long_palette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                    "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                    "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                    "#8A7C64", "#599861")
    palette <- c("black", rep(long_palette, ceiling(max(Muller_df$Identity) / length(long_palette))))
  }
  if(is.na(colour_by)) colour_by <- "Identity"
  ggplot(Muller_df, aes(x=Generation, y=Frequency, group = Group_id, fill = as.factor(get(paste0(colour_by))), colour = as.factor(get(paste0(colour_by))))) + 
    geom_area(size = 0.5) + # add lines to conceal the gaps between areas
    scale_fill_manual(values = palette, name = colour_by) + 
    scale_color_manual(values = palette) + 
    theme(legend.position = ifelse(add_legend, "right", "none")) +
    guides(linetype=FALSE,color=FALSE)
}
