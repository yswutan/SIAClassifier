# SIAClassifier
An R package for molecular subtyping of small intestinal adenocarcinoma

Install package
```
devtools::install_github("yswutan/SIAClassifier")
library(SIAClassifier)
```

How to use
```
## Expr: a dataframe with Gene Expression Profiles data values,
##       samples in columns, genes in rows, rownames corresponding to gene symbols


SIAlabel <- SIAClassifier(Expr = Expr, minPosterior = 0.5, scale=TRUE, log2transfrom=FALSE)
```
