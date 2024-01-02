multigroup_bd <- function (modelList, group, standardize = "scale", standardize.type = "latent.linear", 
          test.type = "III") 
{
  name <- deparse(match.call()$modelList)
  data <- modelList$data
  modelList <- removeData(modelList, formulas = 1)
  intModelList <- lapply(modelList, function(i) {
    rhs2 <- ifelse(class(i) == "lmerMod", paste(paste(paste(all.vars_trans(i)[-1], "*", group), collapse = " + "),  "+ (", findbars(formula(i)),")"), paste(paste(all.vars_trans(i)[-1], "*", group), collapse = " + "))
    i <- update(i, formula(paste(". ~ ", rhs2)))
    return(i)
  })
  newModelList <- lapply(unique(data[, group]), function(i) update(as.psem(modelList), 
                                                                   data = data[data[, group] == i, ]))
  names(newModelList) <- unique(data[, group])
  coefsList <- lapply(newModelList, coefs, standardize, standardize.type, 
                      test.type)
  names(coefsList) <- unique(data[, group])
  coefTable <- coefs(modelList, standardize, standardize.type, 
                     test.type)
  anovaTable <- anova(as.psem(intModelList))[[1]]
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), 
  ]
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Response", 
                                                   "Predictor")]
  global$Predictor <- sub(":", "\\1", sub(group, "\\1", global$Predictor))
  if (nrow(global) == nrow(anovaInts)) 
    newCoefsList <- list(global = coefTable)
  else {
    newCoefsList <- lapply(names(coefsList), function(i) {
      ct <- as.matrix(coefsList[[i]])
      idx <- which(apply(ct[, 1:2], 1, paste, collapse = "___") %in% 
                     apply(global[, 1:2], 1, paste, collapse = "___"))
      ct[idx, ] <- as.matrix(coefTable[idx, ])
      ct <- cbind(ct, ifelse(1:nrow(ct) %in% idx, "c", 
                             ""))
      for (j in 1:nrow(ct)) {
        if (ct[j, ncol(ct)] == "c") {
          model <- modelList[[which(sapply(listFormula(modelList), 
                                           function(x) all.vars.merMod(x)[1] == ct[j, 
                                                                                   "Response"]))]]
          data. <- data[data[, group] == i, ]
          sd.x <- GetSDx(model, modelList, data., standardize)
          sd.x <- sd.x[which(names(sd.x) == ct[j, "Predictor"])]
          sd.y <- GetSDy(model, data., standardize, standardize.type)
          new.coef <- as.numeric(ct[j, "Estimate"]) * 
            (sd.x/sd.y)
          ct[j, "Std.Estimate"] <- ifelse(length(new.coef) > 
                                            0, round(as.numeric(new.coef), 4), "-")
        }
      }
      ct <- as.data.frame(ct)
      ct[is.na(ct)] <- "-"
      names(ct)[(ncol(ct) - 1):ncol(ct)] <- ""
      return(ct)
    })
    names(newCoefsList) <- names(coefsList)
  }
  if (nrow(global) == nrow(anovaInts)) 
    gof <- fisherC(modelList)
  else {
    b <- basisSet(modelList)
    cf <- coefTable[coefTable$Response %in% global$Response & 
                      coefTable$Predictor %in% global$Predictor, ]
    b <- lapply(b, function(i) {
      for (j in 3:length(i)) {
        value <- cf[cf$Response == i[2] & cf$Predictor == 
                      i[j], "Estimate"]
        if (length(value) != 0) 
          i[j] <- paste0("offset(", value, "*", i[j], 
                         ")")
      }
      return(i)
    })
    if (length(b) == 0) 
      b <- NULL
    gof <- fisherC(modelList, basis.set = b)
  }
  ret <- list(name = name, group = group, global = global, 
              anovaInts = anovaInts, group.coefs = newCoefsList, Cstat = gof)
  class(ret) <- "multigroup.psem"
  return(ret)
}
