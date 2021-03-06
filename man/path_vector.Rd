% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MullerPlot.R
\name{path_vector}
\alias{path_vector}
\title{Record a path through all nodes of an adjacency matrix}
\usage{
path_vector(edges)
}
\arguments{
\item{edges}{Dataframe comprising an adjacency matrix, with column names "Parent" and "Identity"}
}
\value{
A vector specifying the path.
}
\description{
Nodes are traversed in the order that they should be stacked in a Muller plot. 
Each node appears exactly twice.
}
\examples{
edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
path_vector(edges1)

}
\author{
Rob Noble, \email{robjohnnoble@gmail.com}
}
