#' Title Gene expression data imputation
#'
#' @param Exp a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols
#' @param G a character containing gene symbol(s)
#'
#' @return a Gene Expression Profiles after imputing NA values
#' @export
naImpute <- function(Exp,G=NULL){
  m <- rowMeans(Exp,na.rm=T)
  mm <- mean(m,na.rm=T)
  m[which(is.na(m))] <- mm
  mna <- rowMeans(Exp)
  w <- which(is.na(mna))
  for(i in w){
    try(Exp[i,which(is.na(Exp[i,]))] <- m[i])
  }
  G <- setdiff(G,rownames(Exp))
  if(!is.null(G)) {
    tmp <- as.data.frame(array(mm,dim=c(length(G),ncol(Exp))))
    rownames(tmp) <- G
    names(tmp) <- names(Exp)
    Exp <- rbind(Exp,tmp)
  }
  Exp
}
#' Title small intestinal adenocarcinoma molecular subtypes prediction
#'
#' @param Expr a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols
#' @param minPosterior minimal posterior probability to classify a sample. By default is 0.5.
#' @param scale scaling would be performed (when scale=TRUE)
#' @param log2transfrom log2transfrom would be performed (when log2transfrom=TRUE)
#'
#' @return a data frame with SIA subtype label and probability for each sample.
#' @export
#'
#' @import randomForest
SIAClassifier <- function(Expr = NULL, minPosterior = 0.5, scale=TRUE, log2transfrom=FALSE){
  rownames(Expr) <- gsub("-", "_", rownames(Expr))
  featureGenes <- rownames(finalModel$importance)
  if(length(intersect(featureGenes, rownames(Expr))) == 0){
    stop('Rownames should be gene symbol in expression profile.')
  }
  if(log2transfrom){
    Expr_log2 <- log2(Expr)
  } else {
    Expr_log2 <- Expr
  }
  if(scale){
    Expr_scale <- as.data.frame(scale(t(Expr_log2)))
  } else {
    Expr_scale <- as.data.frame(t(Expr_log2))
  }
  if(length(intersect(featureGenes, colnames(Expr_scale))) < length(featureGenes)){
    message('Imputing gene expression profile......')
    Expr_impute <- as.data.frame(t(naImpute(as.data.frame(t(Expr_scale)), featureGenes)))
  } else {
    Expr_impute <- Expr_scale
  }
  label <- as.data.frame(predict(finalModel, newdata= Expr_impute, type = "prob"))
  label$nearest <- predict(finalModel, newdata= Expr_impute)
  label$predict <- apply(label[, 1:4], 1, function(x){
    if(max(x) > minPosterior){
      colnames(label)[which.max(x)]
    } else {NA}
  })
  return(label)
}
