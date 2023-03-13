# ExpSigfinder
Scripts for identifying mutational signatures from experimental datasets. This is a fork of xqzou/ExpSigfinder , using some of the underlying functionality.

WIP - no guarantee any of this currently works

## Installation

NB. this is to install the original repo, will update once i work out how to best implement

```{r, eval = FALSE}
# Install the released version from Github
git clone https://github.com/xqzou/ExpSigfinder.git
cd ExpSigfinder
R CMD INSTALL .
```


## Example

NB. this is an example of how one might run the original underlying code from the original repo, will update

```{r, eval = FALSE}
# import the packages 
library("ExpSigfinder")
library(tidyverse)

df <- read.table("example.tsv", sep = "\t", header = T, as.is = T) 
df <- df %>%
  mutate(bg_profile=rowMeans(df[,c("treatment.1","treatment.2","treatment.3")]))
  
# calculate the average number of mutation burden in controls
bg_mean <- sum(df$bg_profile)

# High
Wrap_KOSig(df,"bg_profile",c(("treatment.1","treatment.2","treatment.3"),100, bg_mean,2,"treatment")
sig <- read.table("treatment.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("treatment",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("treatment","_percentage.pdf"))
```
