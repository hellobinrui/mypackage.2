#'Predict method for linear model fits
#'
#'Predicted values based on linear model object.
#'
#'@param object Object of class inheriting from "lm"
#'@param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#'
#'@return predict.lm produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. For type = "terms" this is a matrix with a column per term and may have an attribute "constant".
#'
#'@examples
#'## Predictions
#'x <- rnorm(15)
#'y <- x + rnorm(15)
#'predict(lm(y ~ x))
#'new <- data.frame(x = seq(-3, 3, 0.5))
#'predict(lm(y ~ x), new, se.fit = TRUE)
#'pred.w.plim <- predict(lm(y ~ x), new, interval = "prediction")
#'pred.w.clim <- predict(lm(y ~ x), new, interval = "confidence")
#'
#'@export
#'


predictlm<-function(object, newdata, interval=c("none", "prediction", "confidence"),
                    type = c("response", "terms"), level = 0.95, se.fit = FALSE,
                    weights = 1){
  # set up the values for important variables
  r = object$residuals
  n = length(r)
  p = object$rank
  beta = object$coefficients
  w = object$weights
  rss = sum(ifelse(is.null(w), r^2 , r^2*w))
  df_res = object$df.residual
  res_var = rss/df_res
  pred_var = res_var/weights
  pivt = qr(object)$pivot[seq_len(p)]
  na.act = object$na.action
  xlev = object$xlevels
  type = match.arg(type)
  interval = match.arg(interval)

  # if the object is not a linear model, send out this warning message
  if (!inherits(object, "lm")) {
    warning("It's not a linear model object!")
  }
  # set up the design matrix, if we don't have new data, then the design matrix is the one used in the object
  if (is.na(newdata) || missing(newdata) || is.null(newdata)) {
    X = model.matrix(object)
  } else { # if we have the new data, first check if the model frame is the data class, then let the design matrix be the one of the new data
    new_terms <- delete.response(terms(object))
    mf <- model.frame(new_terms, newdata, na.action = na.act,
                      xlev = xlev)
    data_cl = attr(new_terms , "dataClasses")
    if (!is.null(data_cl)){
      .checkMFClasses(data_cl, mf)
    }
    X <- model.matrix(new_terms, mf)
  }
  # send out a warning message if the rank of object is less than the number of columns of design matrix  X
  if (p < ncol(X) ){ #&& !(missing(newdata) || is.null(newdata)
    warning("prediction from a rank-deficient fit may be misleading")
  }
  predictor <- drop(X[, pivt, drop = FALSE] %*% beta[pivt])
  #if(type == "response" || is.null(type)){
  #return(predictor)
  #}
  if (interval != "none" || se.fit) {
    if(type == "response"){
      if (p > 0) {
        if (missing(newdata) && is.null(w)){
          inv_xr = qr.Q(qr(object))[, seq_len(p), drop = FALSE]
        } else {
          inv_xr = X[, pivt] %*% qr.solve(qr.R(qr(object))[seq_len(p), seq_len(p)])
        }
        index_p = drop(inv_xr^2 %*% rep(res_var, p))
      } else{
        index_p = rep(0, n)
      }
    }
    if (interval == "prediction") {
      if (missing(newdata)){
        warning("new data is missing and prediction is based on current data\n")
      }
      if (missing(newdata) && missing(weights) && !is.null(w)) {
        weights <- w
        warning("assuming prediction variance inversely proportional to weights used for fitting\n")
      }

    }
  }

  if (type == "terms") {
    x = model.matrix(object)
    avx = colMeans(x)
    np = length(attr(x, "assign"))
    if(attr(terms(object), "intercept") > 0){
      assign = list("x" = seq(2, np))
      X <- sweep(X, 2L, avx, check.margin = FALSE)
      const <- sum(avx[pivt] * beta[pivt])
    } else{
      assign = list("x1"= 1, "x2" = seq(2, np))
    }
    nt <- length(assign) # 1 or 2
    predictor <- matrix(ncol = nt, nrow = nrow(X))
    if (interval!="none"){
      index_p <- matrix(ncol = nt, nrow = nrow(X))
      dimnames(index_p) <- list(rownames(X), names(assign))
      inv_r <- qr.solve(qr.R(qr(object))[seq_len(p), seq_len(p)])
    }
    unpivot <- rep(0, ncol(X))
    unpivot[pivt] <- seq_len(p)
    for (i in seq(1, nt, length.out = nt)) {
      id <- assign[[i]]
      predictor[, i]<-X[, id, drop = FALSE] %*% beta[id]
      if(interval !="none"){
        index_p[, i] <- ifelse (any(id > 0),
                                as.matrix(X[, id, drop = FALSE] %*% inv_r[unpivot[id], , drop = FALSE])^2 %*% rep(res_var, p),
                                0)
      }
    }
    dimnames(predictor) <- list(rownames(X), names(assign))
    attr(predictor, "constant") <- ifelse (attr(terms(object), "intercept") > 0, const, 0)
  }
  if(interval!="none"){
    if(interval == "prediction"){
      bound = qt((1 - level)/2, df_res) * sqrt(index_p + pred_var)
    } else if (interval == "confidence"){
      bound = qt((1 - level)/2, df_res) * sqrt(index_p)
    }
    lwr = predictor- bound
    upr = predictor+ bound
    if (type=="terms"){
      return(list(fit = predictor))
    } else if (type =="response"){
      #predictor = cbind(predictor, lwr, upr)
      colnames(predictor) = c("fit")
      return(predictor)
    }
  } else{
    return(predictor)
  }

}

