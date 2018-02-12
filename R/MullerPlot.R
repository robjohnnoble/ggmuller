#' Move to daughter in adjacency matrix
#'
#' Returns the first Identity value in the sorted set of daughters. 
#' When parent has no daughters, returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param parent number or character string specifying whose daughter is to be found
#'
#' @return The daughter's Identity.
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
move_down <- function(edges, parent) {
  if(!(parent %in% edges$Identity) & !(parent %in% edges$Parent)) stop("Invalid parent.")
  daughters <- filter_(edges, ~Parent == parent)$Identity
  if(length(daughters) == 0) return(parent) # if it is not a parent then don't move
  if(is.factor(daughters)) daughters <- levels(daughters)[daughters]
  return(sort(daughters)[1])
}

#' Move to sibling in adjacency matrix
#'
#' Returns the next Identity value among the sorted set of siblings. 
#' When there is no such sibling, returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param identity number or character string specifying whose sibling is to be found
#'
#' @return The sibling's Identity.
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
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
  parent <- filter_(edges, ~Identity == identity)$Parent
  if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
  siblings <- sort(filter_(edges, ~Parent == parent)$Identity)
  siblings <- siblings[which(siblings == identity) + 1]
  if(length(siblings) == 0) return(identity) # if it is the initial genotype then don't move
  if(is.na(siblings)) return(identity) # if it is the initial genotype then don't move
  if(is.factor(siblings)) siblings <- levels(siblings)[siblings]
  return(siblings)
}

#' Move to parent in adjacency matrix
#'
#' Returns the corresponding Parent value. 
#' When there is no parent (i.e. at the top of the tree), returns the input Identity.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#' @param identity number or character string specifying daughter whose parent is to be found
#'
#' @return The Parent value.
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
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
  parent <- filter_(edges, ~Identity == identity)$Parent
  if(length(parent) == 0) return(identity) # if it is the initial genotype then don't move
  if(is.factor(parent)) parent <- levels(parent)[parent]
  return(parent)
}

#' Move to top of adjacency matrix
#'
#' Returns the Parent value of the common ancestor.
#'
#' @param edges Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"
#'
#' @return The Parent that is the common ancestor.
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
  start <- sort(edges$Parent)[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
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
#' path_vector(edges1)
#'
#' @export
#' @import dplyr
path_vector <- function(edges) {
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
    if(n > 2 * dim(edges)[1] + 2) stop("Error: stuck in a loop")
    if(max(table(path) > 2)) stop("Error: adjacency matrix seems to include loops.")
  }
  if(length(path) != 2 * dim(edges)[1] + 2) stop("Error: adjacency matrix seems to be bipartite.")
  return(path)
}

#' Reorder a Muller plot dataframe by a vector
#'
#' @param df Dataframe with column names "Identity", "Parent", and either "Generation" or "Time", in which each Identity appears exactly twice
#' @param vector Vector of Identity values
#'
#' @return The reordered dataframe.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{path_vector}}
#'
#' @examples
#' df <- data.frame(Generation = c(rep(0, 6), rep(1, 6)), 
#'  Identity = rep(1:6,2), Population = c(1, rep(0, 5), 10, rep(1, 5)))
#' df <- rbind(df, df) # duplicate rows
#' require(dplyr)
#' df <- arrange(df, Generation) # put in chronological order
#' edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6) # adjacency matrix
#' path <- path_vector(edges1) # path through the adjacency matrix
#' reorder_by_vector(df, path)
#'
#' @export
#' @import dplyr
reorder_by_vector <- function(df, vector) {
  Generation <- NULL # avoid check() note
  Identity <- NULL # avoid check() note
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(df) && !("Generation" %in% colnames(df))) {
    colnames(df)[colnames(df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # add unique id column to the vector:
  dup <- duplicated(vector)
  vector <- as.data.frame(vector)
  vector$count <- 1:nrow(vector)
  B <- data.frame(Unique_id = apply(vector, 1, function(x)
    if(dup[as.numeric(x["count"])]) paste0(x["vector"], "a") 
    else x["vector"]))
  B <- data.frame(lapply(B, as.character), stringsAsFactors=FALSE)
  vector <- cbind(vector[-ncol(vector)], B)
  
  # useful parameters:
  gens <- unique(df$Generation) # list of unique time points
  n_gens <- length(gens) # number of unique time points
  n_ids <- dim(df)[1] / n_gens # number of unique ids
  
  # add unique id and group id columns to the dataframe:
  dup <- duplicated(df)
  df$count <- 1:nrow(df)
  df$Generation <- as.character(df$Generation)
  df$Identity <- as.character(df$Identity)
  B <- data.frame(t(apply(df, 1, function(x)
    if(dup[as.numeric(x["count"])]) c(paste0(x["Identity"], "a"), paste0(x["Identity"], "a_", x["Generation"])) 
    else c(x["Identity"], paste0(x["Identity"], "_", x["Generation"])))))
  colnames(B) <- c("Group_id", "Unique_id")
  B <- data.frame(lapply(B, as.character), stringsAsFactors=FALSE)
  df$Generation <- as.numeric(df$Generation)
  df <- cbind(df[-ncol(df)], B)
  
  # keep only the unique id column of the vector, repeated n_gens times:
  vector <- rep(vector$Unique_id, n_gens)
  # concatenate the vector entries with the generation values:
  vector <- paste0(vector, "_", rep(gens, each = n_ids))
  
  # reorder the dataframe by the vector:
  df <- df[match(vector, df$Unique_id), ]
  
  # restore original time column name:
  colnames(df)[colnames(df) == "Generation"] <- original_colname
  
  return(df)
}

#' Add rows to a population dataframe to ensure genotype starting points are plotted correctly
#' 
#' 
#' The function 1) identifies when genotypes first have non-zero populations;
#' 2) copies all the rows of data for these time points; 3) modifies the copied rows by decreasing
#' Generation and setting Population of the emerging genotypes to be close to zero;
#' and then 4) adds the modified rows to the dataframe. This ensures that ggplot plots
#' genotypes arising at the correct time points.
#' 
#' By default, the function assumes that each genotype arose half way between the latest time at which 
#' its population is zero and the earliest time at which its population is greater than zero. You can 
#' override this assumption using the start_positions parameter. If start_positions = 0 (respetively 1) 
#' then each genotype is assumed to have arisen at the earliest (respectively latest) time compatible with the data. 
#' Intermediate values are also permitted.
#' 
#'
#' @param pop_df Dataframe with column names "Identity", "Population", and either "Generation" or "Time"
#' @param start_positions Numeric value between 0 and 1 that determines the times at which genotypes are assumed to have arisen (see examples)
#'
#' @return The input Dataframe with additional rows.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#'
#' @examples
#' pop1 <- data.frame(Generation = rep(1:5, each = 4), Identity = rep(1:4, 5), 
#'                    Population = c(1,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1))
#' add_start_points(pop1)
#' 
#' # to see the effect of changing start_positions, compare the Generation columns:
#' add_start_points(pop1, 0)
#' add_start_points(pop1, 1)
#'
#' @export
#' @import dplyr
add_start_points <- function(pop_df, start_positions = 0.5) {
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(pop_df) && !("Generation" %in% colnames(pop_df))) {
    colnames(pop_df)[colnames(pop_df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # set small time interval:
  all_gens_list <- unique(pop_df$Generation)
  delta <- abs(min(1E-2 * min(diff(all_gens_list)), 1E-4 * (max(all_gens_list) - min(all_gens_list))))
  start_positions <- max(start_positions, delta)
  start_positions <- min(start_positions, 1 - delta)
    
  # set small initial population size:
  init_size <- 0
  
  # get reference list of generations at which new genotypes appear (and previous generations):
  min_gen <- min(pop_df$Generation)
  first_gens <- group_by_(pop_df, ~Identity) %>%
    filter_(~max(Population) > 0) %>% 
    summarise_(start_time = ~min(Generation[which(Population > 0)]), 
               previous_time = ~lag(Generation)[min(which(Population > 0))]) %>%
    filter_(~start_time > min_gen) %>%
    ungroup()
  
  # if all genotypes appear at the first time point then don't make any changes:
  if(dim(first_gens)[1] == 0) return(pop_df)
  
  # function to get the generation previous to a specified generation:
  lag_gens <- function(x) {
    ans <- lag(all_gens_list)[which(all_gens_list == x)]
    if(is.na(ans)) return(0)
    return(ans)
  }
  
  # copy all rows for generations at which new genotypes appear:
  gens_list <- unique(first_gens$start_time)
  new_rows <- filter_(pop_df, ~Generation %in% gens_list)
  prev_rows <- filter_(pop_df, ~Generation %in% sapply(gens_list, lag_gens))
  # adjust generations of copied rows:
  new_rows$Generation <- new_rows$Generation - start_positions * (new_rows$Generation - sapply(new_rows$Generation, lag_gens))
  # adjust populations of copied rows:
  new_rows$Population <- (1 - start_positions) * new_rows$Population + start_positions * prev_rows$Population
  # add the copied rows to the dataframe:
  pop_df <- bind_rows(pop_df, new_rows) %>%
    arrange_(~Generation, ~Identity)
  
  # adjust generations in reference list:
  first_gens$Generation <- first_gens$start_time - start_positions * (first_gens$start_time - sapply(first_gens$start_time, lag_gens))
  # set small initial populations in reference list:
  first_gens$Population2 <- init_size

  # replace initial populations in the dataframe with values from reference list:
  pop_df <- merge(pop_df, first_gens, all.x = TRUE)
  pop_df$Population <- ifelse(is.na(pop_df$Population2), pop_df$Population, pop_df$Population2)
  pop_df <- pop_df[, !(names(pop_df) %in% c("Population2", "start_time", "previous_time"))]
  
  # restore original time column name:
  colnames(pop_df)[colnames(pop_df) == "Generation"] <- original_colname

  return(pop_df)
}

#' Create a data frame from which to create a Muller plot
#' 
#'
#' @param edges Dataframe comprising an adjacency matrix, or tree of class "phylo"
#' @param pop_df Dataframe with column names "Identity", "Population", and either "Generation" or "Time"
#' @param cutoff Numeric cutoff; genotypes that never become more abundant than this value are omitted
#' @param start_positions Numeric value between 0 and 1 that determines the times at which genotypes are assumed to have arisen (see examples)
#' @param threshold Depcrecated (use cutoff instead, but note that "threshold" omitted genotypes that never become more abundant than *twice* its value)
#' @param add_zeroes Deprecated (now always TRUE)
#' @param smooth_start_points Deprecated (now always TRUE)
#'
#' @return A dataframe that can be used as input in Muller_plot and Muller_pop_plot.
#'
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{Muller_plot}} \code{\link{Muller_pop_plot}}
#'
#' @examples
#' # by default, all genotypes are included, 
#' # but one can choose to omit genotypes with max frequency < cutoff:
#' Muller_df <- get_Muller_df(example_edges, example_pop_df, cutoff = 0.01)
#'
#' # the genotype names can be arbitrary character strings instead of numbers:
#' example_edges_char <- example_edges
#' example_edges_char$Identity <- paste0("foo", example_edges_char$Identity, "bar")
#' example_edges_char$Parent <- paste0("foo", example_edges_char$Parent, "bar")
#' example_pop_df_char <- example_pop_df
#' example_pop_df_char$Identity <- paste0("foo", example_pop_df_char$Identity, "bar")
#' Muller_df <- get_Muller_df(example_edges_char, example_pop_df_char, cutoff = 0.01)
#'
#' # the genotype names can also be factors (which is the default for strings in imported data):
#' example_edges_char$Identity <- as.factor(example_edges_char$Identity)
#' example_edges_char$Parent <- as.factor(example_edges_char$Parent)
#' example_pop_df_char$Identity <- as.factor(example_pop_df_char$Identity)
#' Muller_df <- get_Muller_df(example_edges_char, example_pop_df_char, cutoff = 0.01)
#'
#' # to see the effect of changing start_positions, compare these two plots:
#' edges1 <- data.frame(Parent = c(1,2,1), Identity = 2:4)
#' pop1 <- data.frame(Time = rep(1:4, each = 4), 
#'                     Identity = rep(1:4, times = 4),
#'                     Population = c(1, 0, 0, 0, 
#'                                    2, 2, 0, 0, 
#'                                    4, 8, 4, 0, 
#'                                    8, 32, 32, 16))
#' df0 <- get_Muller_df(edges1, pop1, start_positions = 0)
#' df1 <- get_Muller_df(edges1, pop1, start_positions = 1)
#' Muller_plot(df0)
#' Muller_plot(df1)
#'
#' @export
#' @import dplyr
#' @importFrom stats na.omit
get_Muller_df <- function(edges, pop_df, cutoff = 0, start_positions = 0.5, threshold = NA, add_zeroes = NA, smooth_start_points = NA) {
  Population <- NULL # avoid check() note
  Generation <- NULL # avoid check() note
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(pop_df) && !("Generation" %in% colnames(pop_df))) {
    colnames(pop_df)[colnames(pop_df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # add missing population values:
  if(dim(pop_df)[1] != length(unique(pop_df$Identity)) * length(unique(pop_df$Generation))) {
    added_rows <- expand.grid(Identity = unique(pop_df$Identity), Generation = unique(pop_df$Generation))
    pop_df <- merge(added_rows, pop_df, all = TRUE)
    pop_df[is.na(pop_df$Population), "Population"] <- 0
    pop_df <- arrange_(pop_df, ~Generation)
    warning("missing population sizes replaced by zeroes")
  }
  
  if (!missing(add_zeroes)) {
    warning("argument add_zeroes is deprecated (it is now always TRUE).", 
    call. = FALSE)
  }
  if (!missing(smooth_start_points)) {
    warning("argument smooth_start_points is deprecated (it is now always TRUE).", 
    call. = FALSE)
  }
  if (!missing(threshold)) {
    warning("argument threshold is deprecated (use cutoff instead, noting that genotypes whose abundance never exceeds the cutoff value are removed, 
            whereas previously genotypes whose abundance never exceeded *twice* the threshold value were removed).", 
            call. = FALSE)
    if (missing(cutoff)) cutoff <- threshold * 2
  }
  
  # check/set column names:
  if(!("Generation" %in% colnames(pop_df)) | !("Identity" %in% colnames(pop_df)) | !("Generation" %in% colnames(pop_df))) 
    stop("colnames(pop_df) must contain Generation (or Time), Identity and Population")
  if("phylo" %in% class(edges)) {
    collapse.singles(edges)
    edges <- edges$edge
  }
  edges <- na.omit(edges) # remove any rows containing NA
  colnames(edges) <- c("Parent", "Identity")
  if(is.factor(edges$Parent)) edges$Parent <- levels(edges$Parent)[edges$Parent]
  if(is.factor(edges$Identity)) edges$Identity <- levels(edges$Identity)[edges$Identity]
  
  # add rows to pop_df to ensure genotype starting points are plotted correctly:
  pop_df <- add_start_points(pop_df, start_positions)
  
  # construct a dataframe with "Age" of each genotype:
  pop_df <- arrange_(pop_df, ~-Population)
  pop_df <- arrange_(pop_df, ~Generation)
  lookup <- group_by_(pop_df, ~Identity) %>% 
    filter_(~Population > 0 | Generation == max(Generation)) %>% 
    slice(1) %>% 
    arrange_(~Generation) %>% 
    ungroup()
  lookup <- mutate(lookup, Age = 1:dim(lookup)[1]) %>% 
    select_(~-c(Generation, Population))
  if(is.factor(lookup$Identity)) lookup$Identity <- levels(lookup$Identity)[lookup$Identity]
  lookup <- select_(lookup, ~c(Identity, Age))
  
  # add semi-frequencies:
  pop_df <- pop_df %>% group_by_(~Generation) %>% 
    mutate(Frequency = (Population / sum(Population)) / 2) %>%
    ungroup()
  pop_df$Population <- pop_df$Population / 2 # because of the duplication
  pop_df$Frequency[is.nan(pop_df$Frequency)] <- 0
  
  # replace each genotype name in adjacency matrix with corresponding Age:
  edges <- filter_(edges, ~Identity %in% lookup$Identity)
  edges <- left_join(edges, lookup, by = "Identity")
  edges <- select_(edges, ~-Identity)
  colnames(edges) <- c("Parent", "Identity")
  edges <- arrange_(edges, ~Identity)
  colnames(lookup)[1] <- "Parent"
  edges <- left_join(edges, lookup, by = "Parent")
  edges$Parent <- edges$Age
  edges <- select_(edges, ~-Age)
  
  # duplicate rows:
  Muller_df <- rbind(pop_df, pop_df)
  Muller_df <- arrange_(Muller_df, ~Generation)
  
  # get the path:
  path <- path_vector(edges)
  path <- rev(path) # apparently, the convention for Muller plots to have earliest-arriving genotypes plotted nearest the top
  
  # replace each Age in the path with corresponding genotype name:
  path <- left_join(data.frame(Age = path), lookup, by = "Age")$Parent
  
  # rearrange the population data according to the path:
  Muller_df <- reorder_by_vector(Muller_df, path)
  
  # optionally remove rare genotypes, and recalculate frequencies:
  if(cutoff > 0) {
    Muller_df <- Muller_df %>% group_by_(~Identity) %>% 
      filter_(~max(Frequency) >= cutoff / 2)
    Muller_df <- Muller_df %>% group_by_(~Generation) %>% 
      mutate(Frequency = Population / sum(Population)) %>% 
      ungroup()
  }
  # the following adjusts for ggplot2 v.2.2.0, which (unlike v.2.1.0) stacks areas in order of their factor levels
  Muller_df$Group_id <- factor(Muller_df$Group_id, levels = rev(
    unlist(as.data.frame(Muller_df %>% filter_(~Generation == max(Generation)) %>% select_(~Group_id)), use.names=FALSE)
    ))
  
  # restore original time column name:
  colnames(Muller_df)[colnames(Muller_df) == "Generation"] <- original_colname
  
  return(Muller_df)
}

#' Draw a Muller plot of frequencies using ggplot2
#'
#' @param Muller_df Dataframe created by get_Muller_df
#' @param colour_by Character containing name of column by which to colour the plot
#' @param palette Either a brewer palette or a vector of colours (if colour_by is categorical)
#' @param add_legend Logical whether to show legend
#' @param xlab Label of x axis
#' @param ylab Label of y axis
#' @param pop_plot Logical for whether this function is being called from Muller_pop_plot (otherwise should be FALSE)
#'
#' @return None
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_Muller_df}} \code{\link{Muller_pop_plot}}
#'
#' @examples
#' # include all genotypes:
#' Muller_df1 <- get_Muller_df(example_edges, example_pop_df)
#' Muller_plot(Muller_df1)
#' # omit genotypes with max frequency < 0.1:
#' Muller_df2 <- get_Muller_df(example_edges, example_pop_df, cutoff = 0.2)
#' Muller_plot(Muller_df2)
#' # colour by a continuous variable:
#' Muller_df1 <- get_Muller_df(example_edges, example_pop_df)
#' Muller_df1$Val <- as.numeric(Muller_df1$Identity)
#' Muller_plot(Muller_df1, colour_by = "Val", add_legend = TRUE)
#' 
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom grDevices col2rgb
Muller_plot <- function(Muller_df, colour_by = "Identity", palette = NA, add_legend = FALSE, xlab = NA, ylab = "Frequency", pop_plot = FALSE) {
  if(!pop_plot & "___special_empty" %in% Muller_df$Group_id) warning("Dataframe is set up for Muller_pop_plot. Use Muller_pop_plot to plot populations rather than frequencies.")
  
  y_factor <- ifelse(pop_plot, "Population", "Frequency")
  if("Time" %in% colnames(Muller_df) && !("Generation" %in% colnames(Muller_df))) x_factor <- "Time"
  else x_factor <- "Generation"
  if(is.na(xlab)) xlab <- x_factor
  direction <- 1
  if(is.na(palette[1])) {
    if(is.numeric(Muller_df[ , colour_by])) {
      palette <- "RdBu"
      direction <- -1
    }
    else {
      long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                      "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                      "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                      "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                      "#5E738F", "#D1A33D")
      palette <- rep(long_palette, ceiling(length(unique(Muller_df$Identity)) / length(long_palette)))
    }
  }
  # test whether palette is a vector of colours; if not then we'll assume it's the name of a predefined palette:
  palette_named <- !min(sapply(palette, function(X) tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)))
  
  gg <- ggplot(Muller_df, aes_string(x = x_factor, y = y_factor, group = "Group_id", fill = colour_by, colour = colour_by)) + 
    geom_area() +
    theme(legend.position = ifelse(add_legend, "right", "none")) +
    guides(linetype = FALSE, color = FALSE) + 
    scale_x_continuous(name = xlab) + 
    scale_y_continuous(name = ylab)
  
  if(is.numeric(Muller_df[ , colour_by])) {
    gg <- gg + 
      scale_fill_distiller(palette = palette, direction = direction, name = colour_by) + 
      scale_color_distiller(palette = palette, direction = direction)
  }
  else {
    if(palette_named) {
      gg <- gg + 
        scale_fill_brewer(palette = palette, name = colour_by) + 
        scale_color_brewer(palette = palette)
    }
    else {
      id_list <- sort(unique(select(Muller_df, colour_by))[[1]]) # list of legend entries, omitting NA
      gg <- gg + 
        scale_fill_manual(values = palette, name = colour_by, breaks = id_list) + 
        scale_color_manual(values = palette)
    }
  }
  return(gg)
}

#' Draw a Muller plot of population sizes using ggplot2
#'
#' This variation on the Muller plot, which shows variation in population size as well as frequency, is also known as a fish plot.
#'
#' @param Muller_df Dataframe created by get_Muller_df
#' @param colour_by Character containing name of column by which to colour the plot
#' @param palette Either a brewer palette or a vector of colours (if colour_by is categorical)
#' @param add_legend Logical whether to show legend
#' @param xlab Label of x axis
#' @param ylab Label of y axis
#'
#' @return None
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_Muller_df}} \code{\link{Muller_plot}}
#'
#' @examples
#' Muller_df <- get_Muller_df(example_edges, example_pop_df)
#' Muller_pop_plot(Muller_df)
#'
#' @export
#' @import dplyr
#' @import ggplot2
Muller_pop_plot <- function(Muller_df, colour_by = "Identity", palette = NA, add_legend = FALSE, xlab = NA, ylab = "Population") {
  
  # add rows for empty space (unless this has been done already):
  if(!"___special_empty" %in% Muller_df$Group_id) Muller_df <- add_empty_pop(Muller_df)
  
  Muller_plot(Muller_df, colour_by = colour_by, palette = palette, add_legend = add_legend, pop_plot = TRUE, xlab = xlab, ylab = ylab)
}

#' Modify a dataframe to enable plotting of populations instead of frequencies
#' 
#' The function adds rows at each time point recording the difference between the total population and its maximum value.
#' Generally there is no need to use this function as Muller_pop_plot calls it automatically.
#'
#' @param Muller_df Dataframe created by get_Muller_df
#'
#' @return A dataframe that can be used as input in Muller_plot.
#' 
#' @author Rob Noble, \email{robjohnnoble@gmail.com}
#' @seealso \code{\link{get_Muller_df}} \code{\link{Muller_pop_plot}}
#'
#' @examples
#' Muller_df <- get_Muller_df(example_edges, example_pop_df)
#' Muller_df2 <- add_empty_pop(Muller_df)
#'
#' @export
#' @import dplyr
#' @import ggplot2
add_empty_pop <- function(Muller_df) {
  Population <- NULL # avoid check() note
  Generation <- NULL # avoid check() note
  . <- NULL # avoid check() note
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(Muller_df) && !("Generation" %in% colnames(Muller_df))) {
    colnames(Muller_df)[colnames(Muller_df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # get maximum total population:
  totals <- Muller_df %>% group_by_(~Generation) %>% 
    summarise_(tot = ~sum(Population)) %>% 
    ungroup
  max_tot <- max(totals$tot)
  
  # avoid warning when Group_id is a factor:
  Muller_df$Group_id <- as.character(Muller_df$Group_id)
  
  # add a new row at start of each Generation group:
  Muller_df <- Muller_df %>%
    group_by_(~Generation) %>%
    summarise_(Identity = ~NA, Population = ~-sum(Population)/2 + 1.1 * max_tot/2) %>%
    mutate(Frequency = NA, 
           Group_id = "___special_empty", 
           Unique_id = paste0("___special_empty_", Generation)) %>% 
    bind_rows(., Muller_df) %>% 
    arrange_(~Generation) %>%
    ungroup()
  
  # add a new row at end of each Generation group:
  Muller_df <- Muller_df %>%
    group_by_(~Generation) %>%
    summarise_(Identity = ~NA, Population = ~first(Population)) %>%
    mutate(Frequency = NA, 
           Group_id = "___special_emptya", Unique_id = paste0("___special_emptya_", Generation)) %>% 
    bind_rows(Muller_df, .) %>% 
    arrange_(~Generation) %>%
    ungroup()
  
  # recalculate frequencies:
  Muller_df <- Muller_df %>% group_by_(~Generation) %>% 
    mutate(Frequency = Population / sum(Population)) %>%
    ungroup()
  
  # the following adjusts for ggplot2 v.2.2.0, which (unlike v.2.1.0) stacks areas in order of their factor levels
  Muller_df$Group_id <- factor(Muller_df$Group_id, levels = rev(
    unlist(as.data.frame(Muller_df %>% filter_(~Generation == max(Generation)) %>% select_(~Group_id)), use.names=FALSE)
  ))
  
  # restore original time column name:
  colnames(Muller_df)[colnames(Muller_df) == "Generation"] <- original_colname
  
  return(Muller_df)
}
