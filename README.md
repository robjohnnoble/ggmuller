ggmuller
========

Create Muller Plots of Evolutionary Dynamics

### Installation

To install & load in R with devtools:

``` r
install.packages("devtools")
library(devtools)
  
install_github("robjohnnoble/ggmuller")
library(ggmuller)
```

### Basic usage

The main functions in `ggmuller` are `get_Muller_df` and `Muller_plot`, which we can run on some data included in the package:

``` r
# get data
Muller_df <- get_Muller_df(example_edges, example_pop_df, threshold = 0.005)

# generate pretty plot
Muller_plot(Muller_df)
```

### How-to guide with examples

For an introduction to how the package works, an overview of its features, and some worked examples, see the blog post on [Visualizing evolutionary dynamics with ggmuller](https://thesefewlines.wordpress.com/2016/08/20/how-to-ggmuller/).

### Version history

0.1.0: First release.

0.1.1: Essential update for compatibility with ggplot2 version 2.2.0.

0.1.2: Better plotting of the points at which genotypes emerge.

0.1.3: Compatibility with dplyr version 0.7.2; add option for smoother plotting of genotype emergence points.

0.2.0: New function Muller_pop_plot shows variation in population size as well as frequency (also known as a fish plot).

### Citation

If you use ggmuller then please ensure your citation includes DOI:10.5281/zenodo.240589.

[![DOI](https://zenodo.org/badge/60275411.svg)](https://zenodo.org/badge/latestdoi/60275411)

