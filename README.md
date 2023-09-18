[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ggmuller)](https://cran.r-project.org/package=ggmuller)

ggmuller
========

Create Muller Plots of Evolutionary Dynamics

### Installation

To install the [CRAN version](https://cran.r-project.org/package=ggmuller):
``` r
install.packages("ggmuller")
library(ggmuller)
```

To install the github version using devtools:

``` r
install.packages("devtools")
library(devtools)
  
install_github("robjohnnoble/ggmuller")
library(ggmuller)
```

### Basic usage

The main functions in `ggmuller` are `get_Muller_df` and `Muller_plot`, which we can run on some data included in the package:

``` r
# reformat data ready for plotting:
Muller_df <- get_Muller_df(example_edges, example_pop_df)

# generate the plot:
Muller_plot(Muller_df)
```

### How-to guides

The best place to start is the [CRAN vignette](https://cran.r-project.org/package=ggmuller/vignettes/ggmuller.html), which includes an overview of features and some worked examples.

An older, slightly different set of instructions can be found in a blog post on [Visualizing evolutionary dynamics with ggmuller](https://thesefewlines.wordpress.com/2016/08/20/how-to-ggmuller/).

### Version history

0.5.6: Avoid loading or importing ape package, to avoid namespace conflict between ape and dplyr.

0.5.5: Optional removal of rare types is now much more efficient for large data sets.

0.5.4: Fix a rare bug; add tip labels to trees.

0.5.3: get_Muller_df now fills all columns when adding missing rows (issue #10); polygon edges now aren't drawn by default. 

0.5.2: Can cope with the special case of only one genotype.

0.5.1: Can use brewer palettes.

0.5: Smoother plots from sparse data; time column can now be called "Time" instead of "Generation"; replace missing population sizes with zeroes.

0.4: Smoother plotting by default when genotypes suddenly grow from zero to large frequencies; new "start_positions" parameter can be used to override this default.

0.3: Bug fixes; in particular, the "threshold" option is replaced by "cutoff" (genotypes whose abundance never exceeds "cutoff" are removed, whereas previously genotypes whose abundance never exceeded *twice* "threshold" were removed).

0.2.2: Bug fixes and deprecation of superfluous options.

0.2.1: Correct mistake in calculating population sizes for Muller_pop_plot.

0.2.0: New function Muller_pop_plot shows variation in population size as well as frequency (also known as a fish plot).

0.1.3: Compatibility with dplyr version 0.7.2; add option for smoother plotting of genotype emergence points.

0.1.2: Better plotting of the points at which genotypes emerge.

0.1.1: Essential update for compatibility with ggplot2 version 2.2.0.

0.1.0: First release.

### Citation

To cite ggmuller in publications please use

    Robert Noble (2019). ggmuller: Create Muller Plots of Evolutionary Dynamics. R package version 0.5.3. doi:10.5281/zenodo.591304 https://CRAN.R-project.org/package=ggmuller

A BibTeX entry for LaTeX users is

    @Manual{,
    title = {ggmuller: Create Muller Plots of Evolutionary Dynamics},
    author = {Robert Noble},
    year = {2019},
    note = {R package version 0.5.3},
    url = {https://CRAN.R-project.org/package=ggmuller},
    doi = 10.5281/zenodo.591304
    }

Please amend the version number as appropriate.
