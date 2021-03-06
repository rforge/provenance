\name{clustree}
\alias{clustree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical elements to plot a dendrogram
}
\description{
	Takes a hclust structure and returns rendered coordinates to draw the tree, according to orientation and type parameters. Also returns leaf labels. Generates segments of tree based on pre-calculated coordinates for leaves and in any orientation, unlike plot.hclust() or package ggdendro.
}
\usage{
clustree<-function(atree,orientation=c(1:4),type=c(1,2),positions=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{atree}{
A tree structure generated by hclust.
}
  \item{orientation}{
branching direction: 1 top-down, 2 right-to-left, 3 bottom-up, 4 left-to-right.
}
  \item{type}{
1 uniform height steps, 2 height calculated by hclust w/base 0, 3 same as 2, w/base = height from hclust.
}
  \item{positions}{
Y-positions in the order of drawing layout (atree$order), i.e. ascending values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Rittner
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotKDE}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
