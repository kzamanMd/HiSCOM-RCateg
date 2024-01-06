mu_prob <- function(cumprob, nobs, ncat1) {
  mat <- matrix(cumprob, nobs, ncat1, TRUE)
  mat <- rbind(mat[, 1], diff(t(mat)), 1 - mat[, ncat1])
  mat <- c(mat)
  return(mat)
}


diag_M <- function(x) {
  x_dims <- length(x)
  mat <- matrix(0, x_dims, x_dims)
  mat[1 + 0:(x_dims - 1) * (x_dims + 1)] <- x
  return(mat)
}

derivmclm <- function(mueta, ncategoriesm1, X_mat) { # nolint
  nobs <- nrow(X_mat) / ncategoriesm1
  ans <- diag_M(rep.int(1, ncategoriesm1))
  ans[seq(2, ncategoriesm1^2, ncategoriesm1 + 1)] <- -1
  ans <- apply(ans, 2, function(x) rep.int(x, nobs))
  mat1 <- matrix(mueta, nobs, ncategoriesm1, TRUE)
  mat1 <- apply(mat1, 2, function(x) rep(x, each = ncategoriesm1))
  mat1 <- ans * mat1
  mat2 <- .rowSums(mat1, nrow(mat1), ncol(mat1), FALSE) *
    X_mat[, -c(1:ncategoriesm1)]
  mat2 <- cbind(mat1, mat2)
  mat2
}



derivmclm_weight <- function(mueta, ncategoriesm1, X_mat) { # nolint
  nobs <- nrow(X_mat) / ncategoriesm1
  ans <- diag_M(rep.int(1, ncategoriesm1))
  ans[seq(2, ncategoriesm1^2, ncategoriesm1 + 1)] <- -1
  ans <- apply(ans, 2, function(x) rep.int(x, nobs))
  mat1 <- matrix(mueta, nobs, ncategoriesm1, TRUE)
  mat1 <- apply(mat1, 2, function(x) rep(x, each = ncategoriesm1))
  mat1 <- ans * mat1
  mat2 <- .rowSums(mat1, nrow(mat1), ncol(mat1), FALSE) * X_mat
  mat2
}


time_dep_dat <- function(F_dat, id, time){
  n <- length(unique(id))
  T <- length(unique(time))
  p <- ncol(as.matrix(F_dat))
  longF <- NULL
  for (i in 1:p){
    newdat <- matrix(0, ncol = T, nrow = n*T)
    for(t in 1:T){##T=3
      newdat[seq(t,n*T, T), t] <- F_dat[seq(t,n*T, T),i]
    }
    longF <- cbind(longF,newdat)
  }
  return(longF)
}


#alpha_hat <- matrix(0, ncat1, ncat1)
#A <- matrix(1,T,T)
#diag(A) <- 0
#FF <- kronecker(A, alpha_hat)

norm    <- function(u) sqrt(sum(u^2)) 
HisCoM_RCateg <- function(y, ID,Time, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(ID))
  y_mat <- as.numeric(Y == Intercept)
  ID <- as.numeric(factor(ID))
  id <- rep(ID, each = ncat - 1)
  Time <- as.numeric(factor(Time))
  time <- rep(Time, each = ncat-1)
  
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  #X.all <- scale(X.all)*sqrt(nobs/(nobs-1))
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(ncase/(ncase-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  
  #############
  #######
  ordindex <- order(id, time)
  y_mat <- y_mat[ordindex]
  X_mat <- X_mat[ordindex, ]
  id <- id[ordindex]
  time <- time[ordindex]
  
  ####
  avec <- as.integer(unlist(lapply(split(id, id), "length")))
  maxclsz <-max(avec)
  maxcl <- maxclsz
  nt <- avec
  #nobs <- sum(nt)
  N <- length(unique(id))
  
  
  
  
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  F_mat1 <- X_mat %*% W2
  indxF <- seq(1,dim(F_mat1)[1],2)
  
  F_mat2 <- time_dep_dat(F_mat1[-indxF,-c(1,2)], ID, Time)
  
  F_mat3 <- apply(F_mat2, 2, function(co) rep(co, each = ncat1))
  F_mat4 <- data.frame(cbind(F_mat1[,c(1,2)], F_mat3))
  colnames(F_mat4) <- paste("X", 1:dim(F_mat4)[2], sep = "")
  F_mat <-  as.matrix(F_mat4, ncol = ncol(F_mat4), dimnames = NULL)
  beta_new <- ginv(t(F_mat) %*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  #wb_new <- W2 %*% beta_new
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat - 1
  
  nbeta <- length(beta_new)
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    #wb_old          <- wb_new
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat1            <- X_mat %*% W2 #### write new way
    F_mat2 <- time_dep_dat(F_mat1[-indxF,-c(1,2)], ID, Time)
    
    F_mat3 <- apply(F_mat2, 2, function(co) rep(co, each = ncat1))
    F_mat4 <- data.frame(cbind(F_mat1[,c(1,2)], F_mat3))
    colnames(F_mat4) <- paste("X", 1:dim(F_mat4)[2], sep = "")
    F_mat <-  as.matrix(F_mat4, ncol = ncol(F_mat4), dimnames = NULL)
    eta              <- drop(F_mat %*% beta_old)   
    
    fitproball       <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob          <- fitproball[-exclude]
    dummy            <- mu.eta(eta)
    pi               <- fitprob
    resids           <- y_mat - pi
    
    T <- length(unique(time))
    #    D_mat <- derivmclm(dummy, ncat1, X_mat)
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      #Dw_mat <- derivmclm_weight(dummy, ncat1, B_mat)
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:N){
        selector         <- id == unique(id)[i] 
        #nn               <- sum(selector)/ncat1
        ans <- diag(ncat1)
        ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
        ans   <- rep(list(ans),T)
        
        ans <- as.matrix(bdiag(ans))
        
        mueta_i1         <-  dummy[selector] 
        mueta_i          <-  apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = T*ncat1))
        dpi_eta          <-  ans*mueta_i
        
        B_mat_i          <-  B_mat[selector,]
        #D_mat            <-  dpi_eta %*% B_mat[selector,]
        
        ####
        Pi_vct           <-  pi[selector]
        V_mat            <-  diag(Pi_vct) - Pi_vct %*% t(Pi_vct) #diag_M(var.pi[selector])  
        V_inv            <-  ginv(V_mat)  #inversematrix
        res_i            <-  resids[selector]
        G_i              <-  t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i              <-  eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i             <-  Z_i - z2[selector]
        Sw_300           <-  (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat           <-  Sw_mat + Sw_300
        Hw_300           <-  (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat           <-  Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j                <-   sqrt(ncase)*w_j/ norm(X_j%*%w_j) 
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat1[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    #diag(t(F_mat)%*% F_mat)
    ###beta_update
    
    F_mat2 <- time_dep_dat(F_mat1[-indxF,-c(1,2)], ID, Time)
    
    F_mat3 <- apply(F_mat2, 2, function(co) rep(co, each = ncat1))
    F_mat4 <- data.frame(cbind(F_mat1[,c(1,2)], F_mat3))
    colnames(F_mat4) <- paste("X", 1:dim(F_mat4)[2], sep = "")
    F_mat <-  as.matrix(F_mat4, ncol = ncol(F_mat4), dimnames = NULL)
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:N){
      selector         <- id == unique(id)[i] 
      ###
      ans <- diag(ncat1)
      ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
      ans   <- rep(list(ans),T)
      ans <- as.matrix(bdiag(ans))
      
      mueta_i1         <-  dummy[selector] 
      mueta_i          <-  apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = T*ncat1))
      dpi_eta          <-  ans*mueta_i
      
      F_mat_i          <-  F_mat[selector,]
      #D_mat            <-  dpi_eta %*% F_mat[selector,]   
      
      
      Pi_vct <-   pi[selector]
      V_mat <-    diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <-    ginv(V_mat)
      res_i <-    resids[selector]
      G_i <-      t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <-      eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <-   (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <-   Sb_mat + Sb_300
      Hb_300 <-   (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <-   Hb_mat + Hb_300
    }
    
    p_beta   <-   rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    beta_new <-   ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new  <-  c(W_new, beta_new)
    #crit <- sum(abs(est_new - est_old))
    
    #wb_new <- W2 %*% beta_new
    #crit <- max(abs(abs(as.numeric(wb_new) - as.numeric(wb_old))))
    crit <- max(abs(as.numeric(est_new) - as.numeric(est_old)))
    
    iter <- iter + 1
    if(iter%%5==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  
  
  ##std calculate##
  var_cov <- ginv(Hb_mat + diag(p_beta)) %*%  Hb_mat %*% ginv(Hb_mat + diag(p_beta))
  std.beta <- sqrt(diag(var_cov))
  ## deviance calculate
  #eta_F <- drop(X_mat %*% W2 %*% beta_coef)
  #fitproball <- mu_prob(linkinv(eta_F), nobs, ncat1)  
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- as.numeric(beta_coef)
  fit$weight_coef <- weight_coef
  #fit$sd <- sdP
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit$std_beta <- std.beta
  fit
}





################
Hiscom_Ord_GEE_CV <- function(y, path, path_var, data, lambda.vec, nfolds){
  n <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  nlambda <- length(lambda.vec)
  #############
  folds = rep(1, n/nfolds)
  for (i in 2:nfolds) {
    folds=c(folds,rep(i,n/nfolds))
  }
  
  ###########
  CV_dat <- NULL
  for (i in 1:nlambda){
    tryCatch({
      cat("NCV", i, "\n")
      lam1 <- lam2 <- lambda.vec[i]
      auc <- c()
      Mse <- c()
      like_Li <- c()
      dev <- c()
      for( k in 1:nfolds){
        tryCatch({
          cat("fold=", k, "\n")
          testindex	<-  which(folds==k,arr.ind=TRUE)
          data_test	= data[testindex,] 
          data_train	= data[-testindex,]  
          ### # Specify Y 
          ytest		= y[testindex]  
          ntest <- length(ytest)
          ytrain		= y[-testindex] 
          ###
          Ridge_model <- HisCoM_RCateg(ytrain, ID = data_train$IID, Time = data_train$Time, path, path_var, 
                                             indx= NULL, data = data_train,
                                             maxiter = 500, lambda1= lam1, lambda2 = lam2,
                                             tol = 0.0001)
          
          w_train <- Ridge_model$W
          beta_train <- Ridge_model$beta_coef
          crit <- Ridge_model$crit
          cat("crit=", crit, "\n")
          ###
          ncase <- length(ytest)
          ID_test <- data_test$IID
          ###
          X.all <- data_test[,match(path_var, colnames(data_test))]
          xnames <- colnames(X.all)
          sex_indx <- which(colnames(X.all) == "sex") 
          X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
          X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
          
          for(i in 1:ncol(X.all)){
            if(i!=sex_indx){
              X_mat_1[,i] <- scale(X_mat_1[,i])*sqrt(ncase/(ncase-1))
            }else{
              X_mat_1[,i] <- X_mat_1[,i]
            }
            
          }
          
          #X_mat_1 <- scale(X_mat_1)*sqrt(ncase/(ncase-1))
          #diag(t(X_mat_1)%*%X_mat_1)
          Intercept1 <- rep.int(seq(ncat-1), length(ID_test))
          X_mat_0 <- model.matrix(~factor(Intercept1)-1 )
          X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
          
          X_mat <- cbind(X_mat_0, X_mat_1)
          
          #####
          eta_F <- drop(X_mat %*% w_train %*% beta_train)
          family <- make.link("logit") 
          linkinv <- family$linkinv
          mu.eta <- family$mu.eta
          nobs <- length(ytest)
          
          
          fitprob <- mu_prob(linkinv(eta_F), nobs, ncat1)  
          
          
          family <- make.link("logit") 
          linkinv <- family$linkinv
          mu.eta <- family$mu.eta
          
          
          fitprob <- mu_prob(linkinv(eta_F), nobs, ncat1)  
          
          
          
          ##
          fitprob1 <- matrix(fitprob, ncol = ncat, byrow=T)
          colnames(fitprob1) <- factor(1:ncat)
          ##3
          pred <- fitprob1
          auc1 <- multiclass.roc(factor(ytest), pred)
          auc[k] <- as.numeric(auc1$auc) 
          
          
          
          Y1 <- rep(ytest, each = ncat)
          Intercept <- rep.int(seq(ncat), length(ytest))
          Y1_mat <- as.numeric(Y1 == Intercept)
          Mse[k] <-(1/length(Y1_mat))*sum((Y1_mat - fitprob)**2)
          like_Li2 <- sum(Y1_mat * log(fitprob))
          dev2 <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 1e-08)) 
          
          tol <- 0.0001
          if(crit < tol){
            dev1 <- dev2
          }else{
            dev1 <- "NA"
          }
          
          if(crit < tol){
            like_Li1 <- like_Li2
          }else{
            like_Li1 <- NA
          }
          cat("dev1=", dev1, "like_Li1 =", like_Li1, "\n")
          cat("dev1=", dev1, "like_Li1 =", like_Li1, "AUC=", auc , "MSE=", Mse , "\n")
          like_Li[k]    <- like_Li1 
          dev[k]        <- dev1
        }, error = function(e){})   
      } ###end K
      
      cv.like_Li <- sum(like_Li)/nfolds
      cv.dev <- sum(dev)/nfolds
      cv.auc <- sum(auc)/nfolds
      cv.MSE <- sum(Mse)/nfolds
      ##
      CV_rslt <- c(lam1, lam2, cv.like_Li, cv.dev, cv.auc, cv.MSE  )
      
      
      cat("CV_Value",  CV_rslt, "\n" )
      CV_dat <- rbind(CV_dat, CV_rslt)
      colnames(CV_dat) <- c("lambda1", "lambda2", "cv.like_Li", "cv.dev", "cv.auc", "cv.MSE")
    }, error = function(e){})
  } ##end i
  return(CV_dat)
}

