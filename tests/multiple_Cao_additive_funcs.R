library("MASS")
library(tidyverse)
library("extraDistr")
library(nleqslv)
library(nloptr)
library(splines)
library(gaussquad)

NonHPP_gen <- function(lambdat, censor, lambda.bar) {
  nn <- 0
  while (nn <= 0) {
    nn <- rpois(1, censor * lambda.bar)
    if (nn > 0) {
      tt <- sort(runif(nn, min = 0, max = censor))
      
      ind <- (sapply(tt, lambdat) /  lambda.bar) > runif(nn)
      tt <- tt[ind]
    } else {
      tt <- numeric(0)
    }
    nn <- length(tt)
  }
  
  tt
}


simhomoPoipro <- function(rate, endTime) {
  num <- rtpois(1, rate * endTime, 0)
  return(sort(endTime * runif(num)))
}

# Simulate homogeneous poisson process
simhomoPoipro_zero <- function(rate, endTime) {
  num <- rpois(1, rate * endTime)
  return(sort(endTime * runif(num)))
}


trans_fun <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s))
  }
  return(res)
}



trans_fun_d <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 1))
  }
  return(res)
}

trans_fun_d1o1 <- function(x, s) {
  if (s == 0) {
    res <- rep(1, length(x))
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-1))
  }
  return(res)
}

trans_fun_d12o1 <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 2))
  }
  return(res)
}

trans_fun_d12o2 <- function(x, s) {
  if (s == 0) {
    res <- rep(1, length(x))
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-2))
  }
  return(res)
}



kerfun <- function(xx){
  pmax((1-xx^2)*0.75,0)
}



simAsytransdata <-
  function(mu,
           mu_bar,
           alpha,
           beta,
           s,
           cen,
           nstep = 20) {
    p <- length(beta)
    
    # 1. Generate censoring times
    cen=runif(1,cen,1.5)
    cen <- min(cen,1)  
    
    
    
    # 2. Generate Z(t) and h(t) as a step function
    
    Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
    
    z <-  2*(pnorm(c(mvrnorm(
      1, rep(0, nstep), Sigmat_z
    ))) - 0.5)
    
    left_time_points <- (0:(nstep - 1)) / nstep
    
    # 3. Generate observation time
    z_fun <- stepfun(left_time_points, c(0, z))

    z_cons<- beta[2] * (mean(z)+rnorm(1,0,1) > 0)

    h_fun <- function(x) {
      beta[1] %*% z_fun(x) + 
        z_cons
    }
    
    # 4. Generate failure time
    if (s == 0) {
      lam_fun <- function(tt)
        exp(alpha(tt) + h_fun(tt))
    } else {
      lam_fun <- function(tt)
        (s * (alpha(tt) + h_fun(tt)) + 1) ^ (1 / s)
    }
    u <- runif(1)
    
    fail_time <- nleqslv(cen/2 , function(ttt)
      legendre.quadrature(lam_fun,
                          lower = 0,
                          upper = ttt,
                          lqrule64) + log(u))$x
    
    
    
    X <- min(fail_time, cen)
    
    # 1. Generate R_ik
    
    obs_times <-
      NonHPP_gen(mu,
                 cen, mu_bar)
    
    if (length(obs_times) == 0)
      obs_times <- cen
    
    covariates_obscov <- matrix(c(z_fun(obs_times),rep((z_cons/beta[2]),length(z_fun(obs_times)))),length(z_fun(obs_times)))
    
    # Results
    return(
      tibble(
        X = X,
        delta = fail_time < cen,
        covariates = covariates_obscov,
        obs_times = obs_times,
        censoring = cen
      ) 
    )
  }


lqrule64 <- legendre.quadrature.rules(64)[[64]]

estproc_additive_Cao_mul <- function(data, n, nknots=NULL, norde=NULL, s=NULL, h, pl = 0) {

  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  p_z <- dim(covariates)[2]
  if (is.null(p_z)) p_z <- 1

  lqrulepoints <- 0.5 * lqrule64$x + 0.5
  outf=function(x){x %o% x}
  outmf=function(x,y){x %o% y}
  
  # A
  
  ## A_indi_inner
  dist_tt <- outer(lqrulepoints, obs_times, "-")
  kerval_tt <- kerfun(dist_tt / h) / h * (dist_tt > 0)
  S0_tt <- outer(lqrulepoints, X, "<=") * kerval_tt
  Zbar_tt <- S0_tt %*% covariates / rowSums(S0_tt)
  Zbar_tt[is.na(Zbar_tt)] <- 0
  cov_Zbar_mat <- expand_grid(cov=covariates,Zbar=Zbar_tt)
  cov_m_Zbar <- cov_Zbar_mat$cov - cov_Zbar_mat$Zbar

  
  A_indi <-
    tibble(
      A_indi_inner = t(mapply(
        outmf,
        array_branch(cov_m_Zbar, 1),
        array_branch(cov_Zbar_mat$cov, 1)
      )) * as.vector(S0_tt) * lqrule64$w,
      subjid = rep(as.integer(id), rep(64, length(id))),
      tid = rep(1:64, length(id))
    ) %>% group_by(subjid) %>% summarise(A = t(colSums(A_indi_inner) * 0.5))
  
  A <- matrix(colSums(A_indi$A)/n,ncol=p_z)
  
  
  # B
  
  dist_XX <- outer(X, obs_times, "-")
  kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
  S0_XX <- outer(X, X, "<=") * kerval_XX
  Zbar_XX <- S0_XX %*% covariates / rowSums(S0_XX)
  Zbar_XX[is.na(Zbar_XX)] <- 0
  B_indi <- tibble(B_indi_inner = (covariates-Zbar_XX)*kerval,subid= as.integer(id)) %>% group_by(subid) %>% summarise(B=t(colSums(B_indi_inner)))
  B <- colSums(B_indi$B)/n
  
  # Point est
  
  beta <- solve(A) %*% B
  
  # Sigma
  
  betaextmat <- do.call(rbind,lapply(beta,function(bb) diag(x=bb,p_z,p_z)))
  A_indi_beta <- A_indi$A %*% betaextmat
  
  Sigma <- matrix(rowSums(apply(B_indi$B-A_indi_beta,1,outf))/n^2,nrow=2)
  
  
  se=sqrt(diag(solve(A) %*% Sigma %*% solve(A)))
  list(est=beta,se=se,A_est=A,B_est=B,Sigma_est=Sigma)
  
}


estproc_ori_dh <- function(data, n, nknots, norder, s, h, initpara=NULL) {
  gammap<- nknots+norder-1
  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)

  knots <- (1:nknots) / (nknots + 1)
  
  bsmat <-
    
    ns(
      X,
      knots = knots,
      # degree = 4,
      intercept = TRUE,
      Boundary.knots = c(0, 1)
    )
  
  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    res <-
      sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)
    
    res
  }
  
  loglik1_1_d <- function(beta, gamma) {
    alphaBeta <- bsmat %*% gamma +covariates %*% beta
    temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
    res <- as.vector(t(temp1) %*% cbind(covariates, bsmat))
    
    res
  }
  
  loglik2_inner_1 <- function(tt, beta, gamma) {
    
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    alpha_tt <-
      
      ns(  tt,
           knots = knots,
           # degree = 4,
           intercept = TRUE,
           Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-
      trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) 
    res <- ifelse(kerval_tt==0,0,res)* kerval_tt
    res <- ifelse(outer(tt, X, "<"),res,0)
    rowSums(res) 
  }
  
  loglik2_inner_1_d <- function(tt, beta, gamma) {
    
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    bsmat_tt <-
      
      ns(
        tt,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_tt <-
      bsmat_tt %*% gamma %>% as.vector()
    temp1 <- trans_fun_d(outer(alpha_tt,c(covariates %*% beta), "+"), s)
    temp1 <- ifelse(outer(tt, X, "<"),temp1,0)
    temp1 <- ifelse(kerval_tt==0,0,temp1) * kerval_tt
    res <- cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
    
    res
  }
  
  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )
    
    res / n
    
  }
  
  loglik_1_d <- function(beta, gamma) {
    res1 <-
      loglik1_1_d(beta, gamma)
    res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
    res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
    (res1 - res2) / n
    
  }
  
  
  f <- function(xx) {
    -loglik_1(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))])
  }

  
  A <- function(beta, gamma) {
    
    # Get Z_bar
    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <-
      
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX, covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates/ rowSums(S0)
    outf=function(x){x %o% x}
    S2 <- S0 %*% ( t(apply(covariates,1,outf))) / rowSums(S0)
    outerprod <-  (S2 - t(apply(S1,1,outf))) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s)
    ) ^ 2 
    outerprod[is.na(outerprod)] <-0
    matrix(colSums( outerprod *  kerval * delta) / n,ncol(data$covariates))
    
  }
  
  B <- function(beta, gamma) {
    
    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <-
      
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX,   covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates / rowSums(S0)
    outf=function(x){x %o% x}
    outerproducinner1 <- (S1 - covariates) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s))
    outerproducinner1[is.na(outerproducinner1)] <-0
    bb1=outerproducinner1*  kerval * delta * ( (X-obs_times)>0 )
    bb1=tibble(a=bb1 ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colSums(a))
    nc=ncol(data$covariates)
    bbvec=NULL
    for(kk in 1:n)
    {
      bbvec=c(bbvec,c(bb1%>% filter(row_number() == kk))$`colSums(a)`)
    }
    bb1=matrix(bbvec,n)
    
    b_inner <- function(tt, beta, gamma) {
      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      bsmat_tt <-
        
        ns(
          tt,
          knots = knots,
          # degree = 4,
          intercept = TRUE,
          Boundary.knots = c(0, 1)
        )
      
      
      
      
      
      alpha_tt <-
        bsmat_tt %*% gamma %>% as.vector()
      inner <- outer(alpha_tt,   covariates %*% beta  , "+")
      S0_tt <- outer(tt, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(tt)) * kerval_tt
      S1_tt <- S0_tt %*% covariates / rowSums(S0_tt)
      res=NULL
      
      for( jj in 1:nc){
        temp1  <-
          outer(tt, X, "<=") * matrix(trans_fun_d(outer(alpha_tt, covariates %*% beta,"+"), s),length(tt))*
          kerval_tt     *outer(as.vector(S1_tt[,jj]),covariates[,jj],"-") 
        temp1[is.na(temp1)] <-0
        r1=tibble(a=t(temp1) ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colMeans(a))
        mr1=NULL
        for(i in 1:length(tt)){
          mr1 =c(mr1,c(r1%>% filter(row_number() == i))$`colMeans(a)`)
        }
        r1=matrix( mr1,n)
        res=cbind(res,r1)
      }
      res
    }
    
    binner <- b_inner (0.5 * lqrule64$x + 0.5, beta, gamma)
    bb2=NULL
    for( jj in 1:nc){
      b1 <- as.vector(0.5 * rowSums(lqrule64$w * binner[,(1+(nc-1)*64):(64*nc)]))
      bb2=cbind(bb2,b1)
    }
    bb=bb1-bb2
    matrix(colSums(t(apply(bb,1,outf))) / n,ncol(data$covariates))
  }
  

  if(is.null(initpara)) initpara <- rep(0, ncol(data$covariates) + gammap)
  if(s==0) {
    estres <- nloptr(
      x0 = initpara,
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  } else {

    ineqmat <- cbind(covariates,bsmat)
    ineqmat <- ineqmat[(X>obs_times)& (abs(X-obs_times) <= h),]
    estres <- nloptr(
      x0 = initpara,
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      eval_g_ineq = function(xx) -ineqmat %*%xx-1/s,
      eval_jac_g_ineq = function(xx) -ineqmat,
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  }
  
  
  return(estres$solution)
}





hcv<-function(simdata,n,K,s){
  
  test_idk <- lapply(split(sample(1:n,n),rep(1:K,n/K)),sort)
  
  nn=11
  hmin <- simdata %>% filter(delta) %>% filter(X>0.75) %>% mutate(a=X-obs_times) %>% filter(a>0) %>% filter(a==min(a)) %>% pull(a)
  hmin<-max(hmin,n^(-0.6))
  hmax <- max(simdata %>% mutate(diff=X-obs_times) %>% filter(diff>0) %>% group_by(id) %>% mutate(diff1=max(diff)) %>% pull(diff1))
  hmax= min(hmax,n^(-0.3)) 
  hn <- exp(seq(log(hmin),log(hmax),length.out=nn+1)[-1])
  
  
  res <- foreach(hh=hn) %do% {
    foreach(test_idx_one = test_idk) %do% {
      foldkpar=estproc_ori_dh(simdata %>% 
                                filter(!(as.numeric(id) %in% test_idx_one)),n*(1-1/K), 3, 3, s, hh) 
      testing <- simdata %>% filter(as.numeric(id) %in% test_idx_one)
      test_logll <- logll_val(foldkpar,testing,n/K, 3, 3, s, hh, pl = 0)
      tibble(test_logll,bd=hh)
    } %>% bind_rows(.id="fold")
  }
  
  res %>% bind_rows() %>% group_by(bd) %>% summarize(cvloss=mean(test_logll))
  
}

#### Calculate loglikelihood


logll_val <- function(para,data, n, nknots, norder, s, h, pl = 0) {
  gammap<- nknots+norder-1
  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  knots <- (1:nknots) / (nknots + 1)
  
  bsmat <-
    
    ns(
      X,
      knots = knots,
      # degree = 4,
      intercept = TRUE,
      Boundary.knots = c(0, 1)
    )
  
  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    res <-
      sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)
    
    res
  }
  
  
  loglik2_inner_1 <- function(tt, beta, gamma) {
    
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    alpha_tt <-
      
      ns(  tt,
           knots = knots,
           # degree = 4,
           intercept = TRUE,
           Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-
      trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) * kerval_tt
    rowSums(outer(tt, X, "<") * res) 
  }
  
  
  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )
    
    res / n
    
  }
  
  f <- function(xx) {
    -loglik_1(xx[1:ncol(simdata$covariates)], xx[-(1:ncol(simdata$covariates))])
  }

  f(para) 
  
}




library(tidyverse)
library(caret)
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores)
