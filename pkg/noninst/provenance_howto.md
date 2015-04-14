---
title: "Package provenance - HowTo"
author: "Martin Rittner"
date: "12/08/2014"
output: pdf_document
---

# 1. Introduction
Package `provenance` provides functions for the visualisation of typical data utilised in sediment provenance analysis in the geosciences, as well as helper functions aiding data analysis. Plotting takes advantage of the `ggplot2` package, and the output plots are ggplot objects, which allows for further modification with `ggplot2`'s functions and easy saving of the graphics.  
The "workhorse functions" of the package are:

`plotKDE()` provides plotting of kernel density estimates (**KDE**s) of geochronological data for single or suites of samples

`plotMDS()` provides plotting of multidimensional scaling (**MDS**) maps for quick visual identification of trends and similar groups within a set of sample data

Working with provenance data using the `provenance` package in R typically involves three main steps:

1. load data files, filter and reformat input data for use
2. plot visualisation of data
3. save output

Of these steps, 2 usually is an iterative "trial-and-error" process to determine best graphical representation of the data, and 3 is, due to `ggplot2`'s functionality, a trivial call to `ggsave()`. The first item, loading the data, is usually the most involved and requires the user to write their own scripts for the purpose, as this can not easily be stadardised over the wide range of possible data file formats and the individual requirements of each user. This step can be further subdivided into:  
> 1a loading data files  
> 1b reformattiong data  
> 1c adding information for visualisation

All basic steps will be presented in workable examples in section 2 of this document. Section 3 illustrates the effects of different settings of the parameters given to the main plotting functions.

#####

# 2. Basic workflow

## Loading data

### Loading data files

### Reformatting and filtering data

### Adding information for visualisation

## Plotting data

## Saving the output

#####


# 3. Detailed description of parameters




This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r}
summary(cars)
```

You can also embed plots, for example:

```{r, echo=FALSE}
plot(cars)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
