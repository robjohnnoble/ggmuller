## ------------------------------------------------------------------------
library(ggmuller)
Muller_df <- get_Muller_df(example_edges, example_pop_df)
Muller_plot(Muller_df)

## ------------------------------------------------------------------------
Muller_df <- get_Muller_df(example_edges, example_pop_df)
Muller_pop_plot(Muller_df)

## ------------------------------------------------------------------------
edges3 <- data.frame(Parent = paste0("clone_", 
 LETTERS[c(rep(1:3, each = 2), 2, 5)]), 
 Identity = paste0("clone_", LETTERS[2:9]))

## ------------------------------------------------------------------------
# a function for generating exponential growth curves:
pop_seq <- function(gens, lambda, start_gen) c(rep(0, start_gen),
                                               exp(lambda * gens[0:(length(gens) - start_gen)]))
lambda <- 0.1 # baseline fitness
gens <- 0:150 # time points
fitnesses <- c(1, 2, 2.2, 2.5, 3, 3.2, 3.5, 3.5, 3.8) # relative fitnesses of genotypes
pop3 <- data.frame(Generation = rep(gens, 9),
 Identity = paste0("clone_", LETTERS[rep(1:9, each = length(gens))]),
 Population = c(1E2 * pop_seq(gens, fitnesses[1]*lambda, 0), 
 pop_seq(gens, fitnesses[2]*lambda, 0), 
 pop_seq(gens, fitnesses[3]*lambda, 10), 
 pop_seq(gens, fitnesses[4]*lambda, 20),
 pop_seq(gens, fitnesses[5]*lambda, 30),
 pop_seq(gens, fitnesses[6]*lambda, 40),
 pop_seq(gens, fitnesses[7]*lambda, 50),
 pop_seq(gens, fitnesses[8]*lambda, 50),
 pop_seq(gens, fitnesses[9]*lambda, 60)),
 Fitness = rep(fitnesses, each = length(gens)))

## ------------------------------------------------------------------------
Muller_df3 <- get_Muller_df(edges3, pop3)

## ------------------------------------------------------------------------
Muller_plot(Muller_df3, add_legend = TRUE, xlab = "Time", ylab = "Proportion")

## ------------------------------------------------------------------------
library(RColorBrewer)
num_cols <- length(unique(Muller_df3$Fitness)) + 1
Muller_df3$Fitness <- as.factor(Muller_df3$Fitness)
Muller_plot(Muller_df3, colour_by = "Fitness", 
 palette = rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(num_cols)), 
 add_legend = TRUE, xlab = "Time", ylab = "Proportion")

## ------------------------------------------------------------------------
library(ggplot2)
my_palette <- c("grey", "red", "magenta", "orange", "yellow", "blue", "darkcyan")
ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
	geom_area(size = 0.5) + # add lines to conceal the gaps between areas
	theme(legend.position = "right") +
	guides(linetype = FALSE, color = FALSE) + 
	scale_y_continuous(labels = 25 * (0:4), name = "Percentage") +
	scale_fill_manual(name = "Identity", values = my_palette) +
	scale_color_manual(values = my_palette)

## ------------------------------------------------------------------------
Muller_df_pop <- add_empty_pop(Muller_df)
ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour = "Identity")) + 
	geom_area(size = 0.5) + # add lines to conceal the gaps between areas
	theme(legend.position = "right") +
	guides(linetype = FALSE, color = FALSE) + 
	scale_fill_manual(name = "Identity", values = my_palette) +
	scale_color_manual(values = my_palette) +
	theme_classic()

## ------------------------------------------------------------------------
tree <- adj_matrix_to_tree(edges3)

## ------------------------------------------------------------------------
library(ape)
tree$tip.label <- 1:length(tree$tip.label) # optional
tree$node.label <- (length(tree$tip.label) + 1):10 # optional
plot(tree, show.node.label = TRUE, show.tip.label = TRUE, tip.color = "red")

