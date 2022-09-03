install.packages("rjson")
install.packages("jsonlite")
library(rjson)
library(jsonlite)

update.packages(ask=FALSE,
                checkBuilt=TRUE,
                repos="https://cloud.r-project.org")
install.packages(c("Rcpp", "caret", "forecast", "ggplot2", "quadprog"), 
                 dependencies=TRUE,
                 repos="https://cloud.r-project.org")



setwd("~/xjj/PORI/oncokb")
data<-jsonlite::stream_in(file("actionableGenes.json"),pagesize = 100)
str(data)
head(data)
