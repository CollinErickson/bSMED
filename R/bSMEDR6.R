if (F) {
  if (!exists('lib.loc')) {lib.loc <- NULL}
  #source("adaptconcept_helpers.R")
  #source('LHS.R')
  #source("random_design.R")
  library(TestFunctions, lib.loc = lib.loc)
  library(cf, lib.loc = lib.loc)
  #library(SMED, lib.loc = lib.loc)
  #library(sFFLHD, lib.loc = lib.loc)
  library(UGP, lib.loc = lib.loc)
  library(magrittr)
  #setOldClass("UGP")
}

#' Class providing object with methods for bSMED
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @importFrom stats optim
#' @keywords data, experiments, adaptive, sequential, simulation, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for running a bSMED experiment.
#' @format \code{\link{R6Class}} object.
#' @examples
#' a <- bSMED$new(D=2,L=1003,func=TestFunctions::gaussian1, obj="func",
#'      n0=0,b=3, nb=5, X0=lhs::maximinLHS(20,2), Xopts=lhs::maximinLHS(500,2))
#' a$run()
#' @field X Design matrix
#' @field Z Responses
#' @field b batch size to select each iteration
#' @field nb Number of batches, includes the first iteration taking X0
#' @field D Dimension of data
#' @field Xopts Available points
#' @field X0 Initial design
#' @field package Which GP package to use in IGP
#' @field stats List of tracked stats
#' @field iteration Which iteration
#' @field mod The GP model from IGP
#' @field obj Character telling what the objective function is, e.g. max or grad
#' @field obj_func Function to evaluate the objective
#' @field obj_alpha A parameter that can be used to have model focus on uncertain areas when model is poor
#' @field scale_obj If the objective is not in [0,1], then this will scale it so it works with this method. Implementation not working now.
#' @field n0 How many points were in initial design, not important since X0 must be given in
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/bSMED}
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the model, adding new data if given, then running optimization again.}
#'   }
bSMED <- R6::R6Class(classname = "bSMED",
  public = list(
    func = NULL, # "function", The function to calculate new values after selected
    func_run_together = NULL, # Should the matrix of values to be run be passed to func as a matrix or by row?, useful if you parallelize your own function or call another program to get actual values
    func_fast = NULL,
    D = NULL, # "numeric",
    L = NULL, # "numeric",
    g = NULL, # "numeric", # g not used but I'll leave it for now
    X = NULL, # "matrix", Z = "numeric", Xnotrun = "matrix",
    #Xnotrun = NULL,
    Z = NULL,
    # s = NULL, # "sFFLHD" an object with $get.batch to get batch of points
    #design = NULL,
    stats = NULL, # "list",
    iteration = NULL, # "numeric",
    obj = NULL, # "character",
    obj_func = NULL, # "function",
    obj_alpha = NULL,
    scale_obj = NULL,
    n0 = NULL, # "numeric"
    #take_until_maxpvar_below = NULL,
    package = NULL, # "character",
    batch.tracker = NULL, # "numeric",
    #force_old = NULL, # "numeric",
    #force_pvar = NULL, # "numeric",
    #useSMEDtheta = NULL, # "logical"
    mod = NULL,
    desirability_func = NULL, # args are mod and XX

    # adding for bSMED
    phat = NULL,
    alpha = NULL,
    use_alpha = NULL,
    alpha_regression_i_values = NULL,
    alpha_regression_p_values = NULL,
    gamma = NULL,
    gammamax = NULL,
    tau = NULL,
    kappa = NULL,
    delta = NULL,
    Xopts=NULL,
    X0=NULL,
    Xopts_removed = NULL,
    pX=NULL,
    pXopts=NULL,
    qX=NULL,
    qXopts=NULL,
    b = NULL,
    p = NULL,
    nb=NULL,
    use_design_region_fix = NULL,
    exchange_algorithm_restarts = NULL,

    parallel = NULL, # Should the new values be calculated in parallel? Not for the model, for getting actual new Z values
    parallel_cores = NULL, # Number of cores used for parallel
    parallel_cluster = NULL, # The object for the cluster currently running

    verbose = NULL, # 0 means print nothing, 1 means print somewhat important things, 2 means print a lot

    initialize = function(D,L=NULL,package="laGP", obj=NULL, n0=0,
                         #force_old=0, force_pvar=0,
                         #useSMEDtheta=F,
                         func, func_run_together=FALSE, func_fast=TRUE,
                         #take_until_maxpvar_below=NULL,
                         #design="sFFLHD",
                         p=NULL, alpha=1, use_alpha=T, scale_obj=NULL,
                         gamma=0, tau=0, kappa=0, use_design_region_fix = FALSE,
                         exchange_algorithm_restarts=1,
                         X0, Xopts, b, nb,
                         estimate.nugget=FALSE, set.nugget=1e-12,
                         parallel=FALSE, parallel_cores="detect",
                         verbose=0,
                         ...) {#browser()
      self$D <- D
      self$L <- L
      self$func <- func
      self$func_run_together <- func_run_together
      self$func_fast <- func_fast
      #self$force_old <- force_old
      #self$force_pvar <- force_pvar
      #self$take_until_maxpvar_below <- take_until_maxpvar_below
      self$verbose <- verbose

      #if (any(length(D)==0, length(L)==0, length(g)==0)) {
      # if (any(length(D)==0, length(L)==0)) {
      if (any(length(D)==0)) {
       message("D must be specified")
      }

      # self$design <- design
      # if (self$design == "sFFLHD") {
      #   #self$s <- sFFLHD::sFFLHD(D=D, L=L, maximin=F)
      # } else if (self$design == "random") {
      #   self$s <- random_design$new(D=D, L=L)
      # } else {
      #   stop("No design 3285729")
      # }
      self$X <- matrix(NA,0,D)
      #self$Xnotrun <- matrix(NA,0,D)
      #if(length(lims)==0) {lims <<- matrix(c(0,1),D,2,byrow=T)}
      #mod$initialize(package = "mlegp")
      # if(is.null(package)) {self$package <- "laGP"}
      # else {self$package <- package}
      self$package <- package
      #self$mod <- UGP$new(package = self$package)
      # Assuming noiseless so I'm setting nugget by default, can be changed by passing in
      self$mod <- UGP::IGP(package = self$package, estimate.nugget=estimate.nugget, set.nugget=set.nugget)
      self$stats <- list(iteration=c(),n=c(),pvar=c(),mse=c(), ppu=c(), minbatch=c(), pamv=c())
      self$iteration <- 1
      self$obj_alpha <- 0.5

      self$Xopts = Xopts
      self$Xopts_removed <- matrix(NA, ncol=self$D, nrow=0)
      self$X0 = X0
      self$b <- b
      # self$qX = NULL,
      # self$qXopts = NULL,


      self$alpha <- alpha
      self$use_alpha <- use_alpha
      self$alpha_regression_i_values <- c()
      self$alpha_regression_p_values <- c()
      self$gamma <- gamma
      self$gammamax <- NA
      self$tau <- tau
      self$kappa <- kappa
      self$delta <- 1e-3 # small positive number, no clue what it should be
      self$nb <- nb
      self$use_design_region_fix <- use_design_region_fix
      self$exchange_algorithm_restarts <- exchange_algorithm_restarts


      self$scale_obj <- scale_obj

      if (is.null(self$scale_obj)) {
       # The objective is scaled to [0,1] if not using alpha and if not told not to.
       # This ensures values are between [0,1], which is needed.
       # Kim et al don't use scaling, this is just an alternative since I don't like their alpha regression.
       self$scale_obj <- !use_alpha
      }

      if (self$scale_obj) { # Scale objective to [0,1]
       #  old version where range is exactly [0,1]
       #  p_scaled_func = function(XX) {#browser()
       #   pred <- self$mod$predict(XX, se=F)
       #   pred2 <- self$mod$predict(matrix(runif(1000*self$D), ncol=self$D), se=F)
       #   predall <- c(pred, pred2)
       #   maxpred <- max(predall)
       #   minpred <- min(predall)
       #   relfuncval <- (pred - minpred) / (maxpred - minpred)
       #   relfuncval
       # }
       p_scaled_func = function(XX) {#print(str(XX));browser()
         pred <- self$mod$predict(XX, se=T)
         XXX <- matrix(runif(1000*self$D), ncol=self$D)
         pred2 <- self$mod$predict(XXX, se=T)
         # predall <- c(pred, pred2)
         maxpred <- max(pred$fit +2*pmax(pred$se, 1e-8), # Adding 1e-8 makes sure all relfuncvals are less than 1
                        pred2$fit+2*pmax(pred2$se, 1e-8))
         minpred <- min(pred$fit, pred2$fit)
         relfuncval <- (pred$fit - minpred) / (maxpred - minpred)
         #print(str(XX))
         #print(summary(relfuncval))
         relfuncval
       }
       self$p <- p_scaled_func
      } else if (self$use_alpha) { # Scale function to [0,1] to avoid problems in calculating q, at least to need to scale below 1
        p_cropped_01 = function(XX) {
          # p values should be between 0 and 1, so fix it here
          pred <- self$mod$predict(XX)
          num_less_0 <- sum(pred < 0)
          num_great_1 <- sum(pred > 1)
          if (num_less_0 > 0) { #browser()
            print(paste(num_less_0, " p values below 0, setting them to 0"))
            pred <- pmax(pred, 0)
          }
          if (num_great_1 > 1) { #browser()
            print(paste(num_great_1, " p values above 1, setting them to 1"))
            pred <- pmin(pred, 1)
          }
          pred
        }
        self$p <- p_cropped_01
      } else { # Otherwise don't scale,
       self$p <- self$mod$predict
      }

      # set objective function to minimize or pick dive area by max
      self$obj <- obj
      if (is.null(self$obj) || self$obj == "mse") { # The default
       #self$obj <- "mse" # Don't want it to be character(0) when I have to check it later
       self$obj_func <- mod$predict.var #function(xx) {apply(xx, 1, mod$predict.var)}
       #function(lims) {
      #   msfunc(mod$predict.var, lims=lims, pow=1, batch=T)
      # }
      } else if (self$obj == "maxerr") {
       self$obj_func <- function(lims) {
         maxgridfunc(self$mod$predict.var, lims=lims, batch=T)
       }
      } else if (self$obj == "grad") {
       self$obj_func <- self$mod$grad_norm#{apply(xx, 1, mod$grad_norm)}
      } else if (self$obj == "func") {
       #self$obj_func <- function(xx) max(1e-16, self$mod$predict(xx))#{apply(xx, 1, mod$grad_norm)}
       #self$obj_func <- function(xx) {pv <- self$mod$predict(xx);ifelse(pv<0,1e-16, pv)}
       self$obj_func <- function(xx) pmax(1e-16, self$mod$predict(xx))
      } else if (self$obj == "pvar") {
       self$obj_func <- function(xx) pmax(1e-16, self$mod$predict.var(xx))#{apply(xx, 1, mod$grad_norm)}
      } else if (self$obj == "gradpvaralpha") {
       self$obj_func <- function(xx) { 1              *      self$mod$grad_norm(xx) +
                                       self$obj_alpha *      max(1e-16, self$mod$predict.se(xx))}
      } else if (self$obj == "nonadapt") {
       # use next batch only #obj_func <<- NULL
      } else if (self$obj == "desirability") {#browser()
       self$obj_func <- function(XX) {list(...)$desirability_func(mod=self$mod, XX=XX)}
       self$desirability_func <- list(...)$desirability_func
      }

      self$n0 <- nrow(X0) #n0
      # Initial with X0 points
      # self$X <- X0 #rbind(self$X, Xnew[1:self$n0, , drop=F])
      # self$Z <- apply(self$X,1,self$func)
      # self$mod$update(Xall=self$X, Zall=self$Z)

      # No longer doing this to initalize data since the points come from a design
      # if (length(self$n0) != 0 && self$n0 > 0) {
      #   Xnew <- matrix(NA, 0, self$D)
      #   while (nrow(Xnew) < self$n0) {
      #     Xnew <- rbind(Xnew, self$s$get.batch())
      #     self$batch.tracker <- rep(self$s$b,self$L)
      #   }
      #   self$X <- rbind(self$X, Xnew[1:self$n0, , drop=F])
      #   self$Z <- c(self$Z, apply(self$X,1,self$func))
      #   self$batch.tracker <- self$batch.tracker[-(1:self$n0)]
      #   #if (nrow(Xnew) > self$n0) {
      #   #  self$Xnotrun <- rbind(self$Xnotrun, Xnew[(self$n0+1):nrow(Xnew), , drop=F])
      #   #}
      #   self$mod$update(Xall=self$X, Zall=self$Z)
      # }

      #if (length(never_dive)==0) {never_dive <<- FALSE}
      #if (length(force_old) == 0) {self$force_old <- 0}
      #if (length(force_pvar) == 0) {self$force_pvar <- 0}
      #self$useSMEDtheta <- if (length(useSMEDtheta)==0) {FALSE} else {useSMEDtheta}

      # Set up parallel stuff
      self$parallel <- parallel
      if (self$parallel) {
       # Use a list to store info about parallel, such as num nodes, cluster, etc
       if (parallel_cores == "detect") {
         self$parallel_cores <- parallel::detectCores()
       } else {
         self$parallel_cores <- parallel_cores
       }

       # For now assume using parallel package

       self$parallel_cluster <- parallel::makeCluster(spec = self$parallel_cores, type = "SOCK")
      }
    },
    run = function(maxit=self$nb - self$iteration + 1, plotlastonly=F, noplot=F) {#browser()
      i <- 1
      while(i <= maxit) {
        #print(paste('Starting iteration', iteration))
        iplotit <- ((i == maxit) | !plotlastonly) & !noplot
        self$run1(plotit=iplotit)
        i <- i + 1
      }

      # If at end, print out results
      if (self$iteration > self$nb) { # At end iteration == nb+1
        self$print_results()
      }
    },
    run1 = function(plotit=TRUE) {#browser()#if(iteration>24)browser()
      if (nrow(self$Xopts) + nrow(self$Xopts_removed) < self$b) {stop("Not enough points left to get a batch #82389, initial design not big enough, b reached")}
      self$update_parameters()
      self$add_data()
      self$update_mod()
      #get_mses()
      #should_dive()
      self$update_stats()
      if (plotit) {
        self$plot1()
      }
      #set_params()
      self$iteration <- self$iteration + 1
    },
    calculate_Z = function(X) {#browser()
      # Used to just be apply(self$X, 1, self$func)
      if (self$parallel && inherits(self$parallel_cluster, "cluster")) {
        # parallel::clusterApply(cl = self$parallal_cluster, x = 1:nrow(X))
        Z <- parallel::parRapply(cl = self$parallel_cluster, x = X, self$func)
      } else if (self$func_run_together) {
        Z <- self$func(X)
      } else {
        Z <- apply(X, 1, self$func)
      }
      if (self$use_alpha) {
        if (any(Z < 0)) {
          warning("A function value is < 0, it is assumed all values are within [0,1]. If it is a near-zero negative number than this should work without major problems.")
        }
        if (any(Z > 1)) {
          warning("A function value is > 1, it is assumed all values are within [0,1]. This may cause major problems.")
        }
      }
      Z
    },
    add_data = function() {#browser()
      # If first time
      if (nrow(self$X) == 0 ) {
        #stop("Shouldn't be here #9238527")
        # self$X <- rbind(self$X, self$s$get.batch())
        self$X <- self$X0
        #self$Z <- c(self$Z,apply(self$X, 1, self$func))
        self$Z <- c(self$Z,self$calculate_Z(X=self$X))
        return()
      }
      #browser()

      # Otherwise use energy to pick points
      if (nrow(self$X) < self$b) {stop("Not enough points to get a full batch")}
      newL <- self$exchange_algorithm() # returns b indices to add




      #;browser()
      # if (self$obj %in% c("nonadapt", "noadapt")) {
      #   Xnew <- self$s$get.batch()
      #   Znew <- apply(Xnew, 1, self$func)
      #   self$X <- rbind(self$X, Xnew)
      #   self$Z <- c(self$Z, Znew)
      #   return()
      # }
      # if (!is.null(self$take_until_maxpvar_below) &&
      #     self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
      #   print(paste("Taking until pvar lower: ", self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
      #   Xnew <- self$s$get.batch()
      #   Znew <- apply(Xnew, 1, self$func)
      #   self$X <- rbind(self$X, Xnew)
      #   self$Z <- c(self$Z, Znew)
      #   return()
      # }

      # Add new points
      # for (iii in 1:5) {
      #   self$Xnotrun <- rbind(self$Xnotrun, self$s$get.batch())
      #   self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, self$L))
      # }
      # newL <- NULL
      #browser()
      # Check if forcing old or pvar
      # if (self$force_old > 0 & self$force_pvar > 0) {
      #   stop("No can force_old and force_pvar")
      # } else if (self$force_old > 0 & self$force_old <= 1) {
      #   rand1 <- runif(1)
      #   if (rand1 < self$force_old) {newL <- 1:self$L}
      # } else if (self$force_old > 1) {
      #   if ((iteration %% as.integer(self$force_old)) == 0) {
      #     newL <- 1:self$L
      #   }
      # } else if (self$force_pvar > 0 & self$force_pvar <= 1) {
      #   rand1 <- runif(1)
      #   if (rand1 < self$force_pvar) {newL <- order(self$mod$predict.var(self$Xnotrun), decreasing=T)[1:self$L]}
      # } else if (self$force_pvar > 1) {
      #   if ((iteration %% as.integer(self$force_pvar)) == 0) {
      #     newL <- order(self$mod$predict.var(self$Xnotrun), decreasing=T)[1:self$L]
      #     #newL <- SMED_selectC(f=mod$predict.var, n=L, X0=X, Xopt=Xnotrun)
      #   }
      # }
      # if nothing forced, run SMED_select
      # if (is.null(newL)) { #browser()
      #   if (F) {# standard min energy
      #     #bestL <- SMED_selectC(f=self$obj_func, n=self$L, X0=self$X, Xopt=self$Xnotrun,
      #     #                      theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      #     Yall <- self$obj_func(rbind(self$X, self$Xnotrun))
      #     Y0 <- Yall[1:nrow(self$X)]
      #     Yopt <- Yall[(nrow(self$X)+1):length(Yall)]
      #     bestL <- SMED_selectYC(n=self$L, X0=self$X, Xopt=self$Xnotrun, Y0=Y0, Yopt=Yopt,
      #                           theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      #     newL <- bestL
      #   } else { # take maximum, update model, requires using se or pvar so adding a point goes to zero
      #     #browser()
      #     gpc <- self$mod$clone()
      #     bestL <- c()
      #     for (ell in 1:self$L) {
      #       #objall <- self$obj_func(rbind(self$X, self$Xnotrun))
      #       objall <- self$desirability_func(gpc, rbind(self$X, self$Xnotrun))
      #       objopt <- objall[(nrow(self$X)+1):length(objall)]
      #       objopt[bestL] <- -Inf # ignore the ones just selected
      #       bestopt <- which.max(objopt)
      #       bestL <- c(bestL, bestopt)
      #       if (ell < self$L) {
      #         Xnewone <- self$Xnotrun[bestopt, , drop=FALSE]
      #         Znewone = gpc$predict(Xnewone)
      #         print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xnotrun)
      #         gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0)
      #       }
      #     }
      #     newL <- bestL#;browser()
      #     #gpc$delete() # This deletes the laGP C side part, don't do it
      #     rm(gpc, objall, objopt, bestopt, bestL, Xnewone, Znewone)#;browser()
      #   }
      # }

      #while(nrow(Xnotrun) < max(L^2, 20)) {
      #if (F) {
      #  objs <- obj_func(Xnotrun)
      #  bestL <- order(objs, decreasing = T)[1:L]
      #} else { # SMED NEW STUFF !!!!
        #browser()
      #  bestL <- SMED_selectC(f=obj_func, n=L, X0=X, Xopt=Xnotrun)
        # cf::cf_func(mod$grad_norm)
        # points(X, col=2, pch=19)
        # text(Xnotrun[,1],Xnotrun[,2])
        # SMED_select(f=obj_func,p=ncol(X),n=8, X0=X, Xopt=Xnotrun)
      #}
      #rand1 <- runif(1)
      #newL <- if (rand1 < force_old) {1:L}
      #        else if (rand1 < force_old + force_pvar) {order(mod$predict.var(Xnotrun), decreasing=T)[1:L]}
      #        else {bestL}#{print(paste('first L',iteration));1:L}

      #Xnew <- self$Xnotrun[newL,]
      #self$Xnotrun <- self$Xnotrun[-newL, , drop=FALSE]
      Xnew <- self$Xopts[newL,]
      self$Xopts <- self$Xopts[-newL, , drop=FALSE]
      self$batch.tracker <- self$batch.tracker[-newL]
      # Znew <- apply(Xnew,1,self$func)
      Znew <- self$calculate_Z(Xnew)
      if (any(duplicated(rbind(self$X,Xnew)))) {browser()}
      self$X <- rbind(self$X,Xnew)
      self$Z <- c(self$Z,Znew)
      #self$update_obj_alpha(Xnew=Xnew, Znew=Znew)
    },
    update_obj_alpha = function(Xnew, Znew) {#browser()
      if (is.null(self$obj_alpha)) return()
      Zlist <- self$mod$predict(Xnew, se.fit=T)
      Zmean <- Zlist$fit
      Zse   <- Zlist$se
      abs.scores <- abs(Znew - Zmean) / Zse
      for (score in abs.scores) {
        if (score < 3) {
          self$obj_alpha <- .5 * self$obj_alpha
        } else {
          self$obj_alpha <- 2  * self$obj_alpha
        }
      }
      print(paste('alpha changed to ', self$obj_alpha))
    },
    update_mod = function() {#browser()
      self$mod$update(Xall=self$X, Zall=self$Z)
    },
    # REMOVED get_mses AND should_dive
    set_params = function() {
    },
    update_stats = function() {#browser()
      # self$stats$ <- c(self$stats$, )
      self$stats$iteration <- c(self$stats$iteration, self$iteration)
      self$stats$n <- c(self$stats$n, nrow(self$X))
      #stats$level <<- c(stats$level, level)
      self$stats$pvar <- c(self$stats$pvar, msfunc(self$mod$predict.var,cbind(rep(0,self$D),rep(1,self$D))))
      self$stats$mse <- c(self$stats$mse, self$mse_func()) #msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
      self$stats$ppu <- c(self$stats$ppu, nrow(self$X) / (nrow(self$X) + nrow(self$Xopts)))
      self$stats$minbatch <- c(self$stats$minbatch, if (length(self$batch.tracker>0)) min(self$batch.tracker) else 0)
      self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
    },
    mse_func = function() {
      if (self$func_fast) {
        msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D)))
      } else {
        NaN
      }
    },
    plot1 = function() {#browser()
      if (self$D == 2) {
        #par(mfrow=c(2,1))
        ln <- 5 # number of lower plots
        split.screen(matrix(
         #c(0,.5,.25,1,  .5,1,.25,1,  0,1/3,0,.25, 1/3,2/3,0,.25, 2/3,1,0,.25),
         c(0,.5,.25,1,  .5,1,.25,1,  0,1/ln,0,.25, 1/ln,2/ln,0,.25, 2/ln,3/ln,0,.25, 3/ln,4/ln,0,.25, 4/ln,1,0,.25),
         ncol=4,byrow=T))
        screen(1)
        #xlim <- lims[1,]
        #ylim <- lims[2,]
        cf::cf_func(self$mod$predict,batchmax=500, pretitle="Predicted Surface ", #pts=X)
              afterplotfunc=function(){points(self$X,pch=19)
                                       points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
              }
        )
        # Plot s2 predictions
        screen(2)
        cf::cf_func(self$mod$predict.var,batchmax=500, pretitle="Predicted Surface ", #pts=X)
               afterplotfunc=function(){points(self$X,pch=19)
                 points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
                 points(self$Xopts, col=2); # add points not selected
               }
        )
        if (self$func_fast) {
          screen(3) # actual squared error plot
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          cf::cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ")
        }
        if (self$iteration >= 2) {
         statsdf <- as.data.frame(self$stats)
         screen(4) # MSE plot
         par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar),
              type='o', log="y", col="white",
              xlab="Iteration", ylab=""
         )
         legend("topright",legend=c("MSE","PVar"),fill=c(1,2))
         points(statsdf$iter, statsdf$mse, type='o', pch=19)
         points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)

         screen(5) # level plot
         #par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         #plot(statsdf$iter, statsdf$minbatch, type='o', pch=19,
         #     xlab="Iteration")#, ylab="Level")
         #legend('bottomright',legend="Batch not run",fill=1)
         cf::cf_func(self$mod$grad_norm, n=20, mainminmax_minmax = F, pretitle="Grad ")

         screen(6) # % of pts used plot
         par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
         legend('bottomleft',legend="% pts",fill=1)
        }
        if (self$func_fast) {
          screen(7) # actual squared error plot
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          cf::cf_func(function(xx){(self$mod$predict(xx) - self$func(xx))^2},
                            n = 20, mainminmax_minmax = F, pretitle="SqErr ")
        }

        close.screen(all = TRUE)
      } else { # D != 2
        par(mfrow=c(2,2))
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        statsdf <- as.data.frame(self$stats)
        #print(ggplot(statsdf, aes(x=iteration, y=mse, col=level)) + geom_line())
        #print(ggplot() +
        #        geom_line(data=statsdf, aes(x=iteration, y=mse, col="red")) +
        #        geom_line(data=statsdf, aes(x=iteration, y=pvar, col="blue"))
        #)
        if (self$iteration >= 2) {
          # 1 mse plot
          plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar),
              type='o', log="y", col="white",
              xlab="Iteration", ylab=""
          )
          legend("topright",legend=c("MSE","PVar"),fill=c(1,2))
          points(statsdf$iter, statsdf$mse, type='o', pch=19)
          points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
          # 2 level plot
          #plot(statsdf$iter, statsdf$level, type='o', pch=19)
          #legend('topleft',legend="Level",fill=1)
          Xplot <- matrix(runif(self$D*50), ncol=self$D)
          if (self$func_fast) {
            Zplot.pred <- self$mod$predict(Xplot)
            Zplot.act <- apply(Xplot,1, self$func)
          } else {
            Zplot.pred <- c()
            Zplot.act <- c()
          }
          Zplot.se <- self$mod$predict.se(Xplot)
          Zused.pred <- self$mod$predict(self$X)
          plot(NULL, xlim=c(min(self$Z, Zplot.act), max(self$Z, Zplot.act)),
              ylim=c(min(Zused.pred, Zplot.pred), max(Zused.pred, Zplot.pred)))
          abline(a = 0, b = 1)
          if (self$func_fast) {points(Zplot.act, Zplot.pred, xlab="Z", ylab="Predicted")}
          points(self$Z, Zused.pred, col=2)
          # 3 % pts used plot
          plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
          legend('bottomleft',legend="% pts",fill=1)
          # 4 grad vs pvar
          Xplot <- matrix(runif(self$D*100), ncol=self$D)
          Xplot_grad <- pmax(1e-8, self$mod$grad_norm(Xplot))#;browser()
          Xplot_se <- pmax(1e-8, self$mod$predict.se(Xplot))
          #if (any(Xplot_se <= 0)) {browser()}
          #if (any(Xplot_grad < 0)) {browser()}
          plot(Xplot_se, Xplot_grad, pch=19, xlab='SE', ylab='Grad', log='xy')
        }
      }
    },
    q = function(X, p=self$p, alpha=self$alpha, gamma=self$gamma) {
      warning("I never use this, maybe should remove, in q #0923848")
      (1 - alpha * p(X)) ^ gamma
    },
    q_from_p = function(p, alpha=self$alpha, gamma=self$gamma) {
      omap <- 1 - alpha * p

      # This section shouldn't be used since it is taken care of
      #  when using alpha regression or scale_obj, only might be
      #  used if using obj as given, which is probably a bad idea.
      # p values should be between 0 and 1, so fix it here
      num_less_0 <- sum(omap < 0)
      num_great_1 <- sum(omap > 1)
      if (num_less_0 > 0) {browser()
        print("Some p values below 0, setting them to 0")
        omap <- pmax(omap, 0)
      }
      if (num_great_1 > 1) { browser()
        print("Some p values above 1, setting them to 1")
        omap <- pmin(omap, 1)
      }

      # (1 - alpha * p) ^ gamma
      omap ^ gamma
    },
    E = function(X1=self$X, X2=NULL) {browser()
      print("Why are you using E? #532355")
      n <- nrow(X1)
      b <- nrow(X2)
      X <- rbind(X1, X2)
      Esum <- 0
      for (i in 1:(n+b-1)) {
        for (j in (i+1):n+b) {
          self$q(X[i,]) * self$q(X[j,]) / sqrt(sum((X[i,] - X[j,])^2))
        }
      }
    },
    Ebn = function(X1=self$X, X2=NULL, q1=NULL, q2=NULL) {#browser()
      n <- nrow(X1)
      b <- nrow(X2)
      X <- rbind(X1, X2) # X2 must be below X1
      qq <- c(q1, q2) # q is quit function
      if (is.null(qq)) {
        print('Calculating q values in Ebn, pass q for faster #5238528')
        qq <- self$q(X)
      }
      Esum <- 0
      for (i in (n+1):(n+b)) {
        for (j in 1:n) {
          #self$q(X[i,]) * self$q(X[j,]) / sqrt(sum((X[i,] - X[j,])^2))
          Esum <- Esum + qq[i] * qq[j] / sqrt(sum((X[i,] - X[j,])^2))
        }
      }
      rm(i)
      for (i in (n+1):(n+b-1)) {
        for (j in (i+1):(n+b)) {
          #self$q(X[i,]) * self$q(X[j,]) / sqrt(sum((X[i,] - X[j,])^2))
          Esum <- Esum + qq[i] * qq[j] / sqrt(sum((X[i,] - X[j,])^2))
        }
      }
      Esum
    },
    Ebn_speedup_lasttwoterms = function(X1=self$X, X2_without_Xl=NULL, Xr, qr, q1=NULL, q2_without_ql=NULL) {#browser()
      n <- nrow(X1)
      b <- nrow(X2_without_Xl)
      X <- rbind(X1, X2_without_Xl) # X2 must be below X1
      qq <- c(q1, q2_without_ql) # q is quit function
      if (is.null(qq)) {
        print('Calculating q values in Ebn, pass q for faster #5238528')
        qq <- self$q(X)
      }
      # Sum over points in design
      Esum3 <- 0
      for (k in 1:n) {
        Esum3 <- Esum3 + qr * qq[k] / sqrt(sum((Xr - X[k,])^2))
      }
      # Sum over candidates temporarily in design
      Esum4 <- 0
      for (k in (n+1):(n+b)) {
        Esum4 <- Esum4 + qr * qq[k] / sqrt(sum((Xr - X[k,])^2))
      }
      # Return both
      Esum3 + Esum4
    },
    exchange_algorithm = function(X=self$X, Xopts=self$Xopts, qX=self$qX, qXopts=self$qXopts, b=self$b, restarts=self$exchange_algorithm_restarts) {#browser()
      best_b_inds <- NULL
      best_Ebn <- Inf
      for (rr in 1:self$exchange_algorithm_restarts) {
        #print(system.time({
          exchange_out <- self$exchange_algorithm_once(X=X, Xopts=Xopts, qX=qX, qXopts=qXopts, b=b)
        #}))
        if (exchange_out$Ebn < best_Ebn) {
          best_b_inds <- exchange_out$b_inds
          best_Ebn <- exchange_out$Ebn
        }
      }
      best_b_inds
    },
    exchange_algorithm_once = function(X=self$X, Xopts=self$Xopts, qX=self$qX, qXopts=self$qXopts, b=self$b) {#browser()
      if (nrow(Xopts) < b) {stop("Not enough points left to select #9813244")}
      if (nrow(Xopts) == b) {
        print("Only b options left, taking those #221771")
        return(1:nrow(Xopts))
      }
      nXopts <- nrow(Xopts)
      b_inds <- sample(1:nXopts, b, replace=FALSE)
      Xb <- Xopts[b_inds, ]
      Ec <- self$Ebn(X1=X, X2=Xb, q1=qX, q2=qXopts[b_inds])

      # The speedup on p16 lets you not recalculate full energy each time
      use_speedup <- TRUE
      # I found the calculations are indeed exact same with and without, so this should be true
      # From total time 2.1 sec to 1.42 on 2D test, so not huge, but could be big reduction for this function
      # Just this function: slow .33 .34 .2  .13
      #                     fast .09 .06 .06 .03, so at least 3x speedup

      # Loop over indices in selected candidate points to replace them if any other candidate lowers energy
      for (l in 1:b) {
        Ebn_star <- Ec
        Ebn_star_index <- NA
        if (use_speedup) {
          Ebnl_twoterms <- self$Ebn_speedup_lasttwoterms(X1=X, X2_without_Xl=Xopts[b_inds[-l],], Xr=Xopts[b_inds[l],], qr=qXopts[b_inds[l]], q1=qX, q2_without_ql=qXopts[b_inds[-l]])
        }
        for (r in setdiff(1:nXopts, b_inds)) {
          # Calculate Ebn where row r places l
          # if (any(is.nan(self$Ebn(X, Xopts[c(b_inds[-l],r),], qX, qXopts[c(b_inds[-l],r)])))) {browser()}
          if (use_speedup) {
            Ebnr_twoterms <- self$Ebn_speedup_lasttwoterms(X1=X, X2_without_Xl=Xopts[b_inds[-l],], Xr=Xopts[r,], qr=qXopts[r], q1=qX, q2_without_ql=qXopts[b_inds[-l]])
            Ebnr <- Ec - Ebnl_twoterms + Ebnr_twoterms
          } else {
            Ebnr <- self$Ebn(X, Xopts[c(b_inds[-l],r),], qX, qXopts[c(b_inds[-l],r)])
          }
          # print(paste("Diff with speedup it", Ebnr_old - Ebnr))
          # If it's the min, update it
          # if (inherits(try(if(Ebnr < Ebn_star){}), "try-error")) {browser()}
          if (Ebnr < Ebn_star) {
            Ebn_star <- Ebnr
            Ebn_star_index <- r
          }
        }
        # If energy is reduced, replace it
        if (Ebn_star < Ec) {
          # Replace index
          b_inds[l] <- Ebn_star_index
          # Update Ec for next round
          Ec <- Ebn_star
        }
      }
      # b_inds now return Ec so it can be run multiple times and take best
      list(Ebn=Ec, b_inds=b_inds)
    },
    update_parameters = function() {
      if (self$iteration == 1) {return()}
      #browser()

      # First get p since it is needed to calculate alpha (might not actually use)
      # pAll <- self$p(rbind(self$X, self$Xopts))
      # self$pX <- pAll[1:nrow(self$X)]
      # self$pXopts <- pAll[-(1:nrow(self$X))]
      # Starting over below to do it right, adding Xoptsremoved to better estimate p
      pAll <- self$p(rbind(self$X, self$Xopts, self$Xoptsremoved))
      self$pX <- pAll[1:nrow(self$X)]
      self$pXopts <- pAll[(1+nrow(self$X)):(nrow(self$X)+nrow(self$Xopts))]

      # update params
      self$kappa <- (self$iteration - 1) / self$nb # p15 right column
      if (self$kappa > 1) {self$kappa <- 1} # Shouldn't be bigger than 1
      # TODO: update self$alpha with model
      self$alpha <- if (self$iteration <=1 || !self$use_alpha) {
                      1
                    } else {
                      self$update_alpha_regression(max_pAll=max(pAll))
                    }
      # print(1 - self$alpha * max(self$p(self$X)))
      #browser()
      self$gammamax <- log(self$delta) / log(1 - self$alpha * max(self$p(self$X))) # p14
      self$gamma <- self$kappa * self$gammamax
      if (self$gammamax==0 || self$gamma==0) {
        warning("gamma and/or gammamax equal zero, so it's not working right. This is very bad and means it isn't doing the energy properly.")
      }

      # Update design region
      if (TRUE & nrow(self$Xopts_removed)>0) { # This add the points back in for consideration
        # don't think it's in paper, but it makes sense since the model changes
        self$Xopts <- rbind(self$Xopts, self$Xopts_removed)
        self$Xopts_removed <- matrix(NA, 0, self$D)
      }
      mod.sample <- self$p(matrix(runif(1e4*self$D), ncol=self$D))
      if (self$use_design_region_fix) {
        # This is an alternative to the paper. It uses the quantile of the average of .05 and .95 instead of the average of the 5th and 95th quantile.
        xi_tau <- (1-self$kappa) * 0.05 + self$kappa * 0.95
        self$tau <- unname(quantile(mod.sample, xi_tau))
      } else {
        xi05 <- unname(quantile(mod.sample, 0.05))
        xi95 <- unname(quantile(mod.sample, 0.95))
        self$tau <- (1-self$kappa) * xi05 + self$kappa * xi95 # p15
      }
      pXopts <- self$p(self$Xopts)
      Xopts_inds_to_remove <- (pXopts < self$tau)
      if (sum(!Xopts_inds_to_remove) < self$b) {#browser()
        print("Tau takes too many, leaving best b")
        Xopts_inds_to_remove <- (pXopts < sort(pXopts, decreasing = TRUE)[self$b])
      }
      if (self$verbose >= 1) {print(paste("Removed", sum(Xopts_inds_to_remove), "points"))}
      self$Xopts_removed <- rbind(self$Xopts_removed, self$Xopts[Xopts_inds_to_remove, ])
      self$Xopts <- self$Xopts[!Xopts_inds_to_remove, ]


      # update p and q values
      pAll <- self$p(rbind(self$X, self$Xopts))
      self$pX <- pAll[1:nrow(self$X)]
      self$pXopts <- pAll[-(1:nrow(self$X))]
      self$qX <- self$q_from_p(self$pX)
      self$qXopts <- self$q_from_p(self$pXopts)
      if (any(is.nan(self$qXopts))) {browser()}
      self$q_from_p(self$pXopts)
      12
    },
    update_p_and_q_values = function() {
      stop("I think this is never used! #8293487372")
      pAll <- self$pgroup(rbind(self$X, self$Xopts))
      self$pX <- pAll[1:nrow(self$X),]
      self$pXopts <- pAll[- 1:nrow(self$X), ]
      self$qX <- self$q_from_p(self$pX)
      self$qXopts <- self$q_from_p(self$pXopts)
    },
    update_alpha_regression = function(max_pAll) {#browser()
      # alpha is used to scale phat so max of alpha*p is at one
      # But if I'm using a scaled relative value, then it doesn't matter
      # So using this with relative values gives issues in log(1-alpha*maxp)
      nb <- self$iteration - 1 # First was initial points, haven't done current yet
      # print("Can't I use initial points?") I think I do now, see if statement below when running iter 2
      # i_values0 <- 1:nb
      # pi_values0 <- sapply(1:nb, function(ii) {max(self$pX[1:(self$n0 + ii * self$b)])})
      # pnb <- pi_values0[nb]


      # p values for alpha regression are calculated after new data is added
      # and model is updated
      # Add value for most recently completed iteration
      self$alpha_regression_i_values <- c(self$alpha_regression_i_values, self$iteration - 1)
      self$alpha_regression_p_values <- c(self$alpha_regression_p_values, max_pAll)
        if (self$iteration == 2) { # if running iter 2, we need to initialize
          self$alpha_regression_i_values <- c(0, self$alpha_regression_i_values)
          self$alpha_regression_p_values <- c(self$alpha_regression_p_values[1]/2, self$alpha_regression_p_values)
        }

      i_values0 <- self$alpha_regression_i_values
      pi_values0 <- self$alpha_regression_p_values
      pnb <- pi_values0[nb+1]

      pnb_bar <- (nb*pnb + 1) / (nb + 1)
      # i_values <- c(i_values0, 0, self$nb)
      # pi_values <- c(pi_values0, pi_values0[1]/2, pnb_bar)
      i_values <- c(i_values0, self$nb)
      pi_values <- c(pi_values0, pnb_bar)

      pgnb_lb <- (pnb + pnb_bar) / 2

      pgnb <- pgnb_lb
      beta0 <- 0
      beta1 <- 0

      pi_hat <- function(i, pgnbhat, beta0hat, beta1hat) {
        pgnbhat * 1/(1 + exp(-(beta0hat + beta1hat*i)))
      }
      optim_fn <- function(par_in) {
        pgnbhat <- par_in[1]
        beta0hat <- par_in[2]
        beta1hat <- par_in[3]
        mean((pi_values - pi_hat(i_values, pgnbhat, beta0hat, beta1hat))^2)
      }

      optim_out <- optim(par=c(pgnb_lb, 0, 0), fn=optim_fn)
      pgnb <- optim_out$par[1]
      # These will plot the data and show the fitted line
      # plot(i_values, pi_values)
      # curve(pgnb * 1/(1 + exp(-optim_out$par[2] - x*optim_out$par[3])), add=T, col=2)

      # If it is below lower bound, set it to be bound
      # The lb is a weighted average of the current pi (pnb) and pnb_bar,
      #  which is between pnb and 1, so lb should be higher than pnb.
      #  But it won't be if pnb is above 1, which can happen if the model is bad?
      #  Should I force p predictions to be between [0,1]???
      if (pnb>1) {
        browser(text = "I think is when there will be trouble, model predicts bigger than 1 but it shouldn't be allowed")
      }
      if (pgnb < pgnb_lb) {
        if (self$verbose >= 1) {
          print("pgnb is lower than lb, setting to be lb to avoid big problems, this is what should be done, trust me (and Kim et al)")
        }
        pgnb <- pgnb_lb
      }
      alpha <- 1/pgnb
      if (self$verbose >=1) {print(paste("new alpha is", alpha))}
      alpha
    },
    print_results = function(give_best_prediction=TRUE) { #browser()
      best_index <- which.max(self$Z)
      bestZ <- self$Z[best_index]
      bestX <- self$X[best_index, ]
      # print(paste0("Best design point is ", paste(bestX, collapse = ''),
      #              " with objective value ", bestZ))
      cat("Best design point is\n\t\t", signif(bestX, 3),
                   "\n\twith objective value\n\t\t", bestZ, '\n')

      if (give_best_prediction) {
        XX <- matrix(runif(1e4*self$D), ncol=self$D)
        ZZ <- self$mod$predict(XX)
        best_index2 <- which.max(ZZ)
        optim(XX[best_index2, ], fn=self$mod$predict)
        optim.out <- optim(XX[best_index2, ], fn=self$mod$predict, lower = rep(0, self$D), upper = rep(1, self$D), method = "L-BFGS-B", control = list(fnscale=-1))
        optim_par <- optim.out$par
        optim_val <- optim.out$val
        cat("Best predicted point over domain is \n\t\t", signif(optim_par, 3),
                                      "\n\twith objective value\n\t\t", optim_val, '\n')
      }
    },
    delete = function() {#browser()
      self$mod$delete()
      if (self$parallel) {
        print("Deleting cluster")
        parallel::stopCluster(cl = self$parallel_cluster)
      }
    }
  )
)

if (F) {
  library(sFFLHD)
  library(UGP)
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(cf)
  require(TestFunctions)
  source('LHS.R')
  library(SMED)

  #gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.01)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  a$run(2)


  #sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}; cf_func(sinumoid)
  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,g=3,func=sinumoid,  obj="grad")
  a$run(10, plotlastonly = T)

  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,g=3,func=RFF_get(), obj="grad")
  a$run(4, plotlastonly = T)

  # higher dim
  a <- adapt.concept2.sFFLHD.R6(D=3,L=8,g=3,func=gaussian1)
  a$run(3)

  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,n0=8,func=banana, obj="grad", force_pvar=.2)
  a$run(1)
  a$run(20, plotl=T)

  # grad cont
  cf::cf_func(a$mod$grad_norm)

  # test run times
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  system.time(a$run(20,plotlastonly = T))
  l <- lineprof::lineprof(a$run(1))
  lineprof::shine(l)

  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="gradpvaralpha", n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD')
  a$run(5)

  # Test desirability function
  des_func <- function(mod, XX) {
    pred <- mod$predict(XX, se=F)
    pred2 <- mod$predict(matrix(runif(1000*2), ncol=2), se=F)
    predall <- c(pred, pred2)
    maxpred <- max(predall)
    minpred <- min(predall)
    des <- (pred - minpred) / (maxpred - minpred)
    des
  }
  des_funcse <- function(mod, XX) {
    pred <- mod$predict(XX, se=T)
    pred2 <- mod$predict(matrix(runif(1000*2), ncol=2), se=T)
    predall <- c(pred$fit, pred2$fit)
    maxpred <- max(predall)
    minpred <- min(predall)
    relfuncval <- (pred$fit - minpred) / (maxpred - minpred)
    des <- 1 + 100 * relfuncval
    des * pred$se
  }
  des_func14 <- function(mod, XX) {
    pred <- mod$predict(XX, se=T)
    des <- apply(XX, 1, function(yy) {if (yy[1] < .5) 4 else 1})
    des * pred$se
  }
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD')
  a$run(5)
  cf::cf(function(x) des_funcse(a$mod, x), batchmax=1e3, pts=a$X)
}
if (F) {
  a <- bSMED$new(D=2,L=1003,func=TestFunctions::gaussian1, obj="grad", n0=0,
                 b=3, nb=4, X0=lhs::maximinLHS(20, 2), Xopts=lhs::maximinLHS(10, 2))
  a$run(2)

}
if (F) {
  quad_peaks <- function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}
  quad_peaks_slant <- TestFunctions::add_linear_terms(function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}, coeffs = c(.02,.01))
  cf::cf(quad_peaks)
  cf::cf(quad_peaks_slant)
  a <- bSMED$new(D=2,L=1003,func=quad_peaks_slant, obj="func", n0=0,b=3, nb=5, X0=rbind(c(.5,.5),lhs::maximinLHS(20,2)), Xopts=lhs::maximinLHS(500,2), use_alpha=T, package="laGP", parallel=F, use_design_region_fix=T, func_fast=T);a$run()
}
