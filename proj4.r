## Aidan Garrity (s1997567), Eleni Michaelidou (s2226022)

## Address of Githup repo: https://github.com/aidangarrity/StatProgramming4.git

## Proportion of the work was undertaken by each team member: 
## Aidan Garrity: Optimization step, error handling, testing
## Eleni Michaelidou: Hessian approximation, robustness checks, testing


## ------------------------- Overview of the code ------------------------------

## Code to create a function to optimize a given objective function w.r.t. some 
## parameters based on the Newton's method. The basic idea of the Newton's method
## is that we start with a guess of the optimization parameter values, and then,
## we evaluate the function and its first (gradient vector) and second derivative
## (Hessian matrix) w.r.t. the parameters at the guess. We use these to obtain
## a quadratic approximation to the objective function. We optimize this approximation,
## and then use it as the parameters for the next guess and approximation for the 
## objective function.
## That process is repeated to convergence. More specifically, convergence is
## achieved when the gradient is approximately zero.

## -----------------------------------------------------------------------------

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,
                 fscale=1,maxit=100,max.half=20,eps=1e-6) {
## Function implementing Newton's method.
## Inputs: theta: a vector of initial values for the optimization parameters.
##                the first argument of func, grad, and hess.
##         func: the objective function to minimize.
##         grad: the gradient function.
##         hess(NULL by default): the Hessian matrix function (if not supplied
##                                then an approximation to the Hessian by finite
##                                differencing of the gradient vector is obtained).
##         toll(1e-8 by default): the convergence tolerance.
##         fscale(1 by default): a rough estimate of the magnitude of func near the optimum.
##         maxit(100 by default): the maximum number of Newton iterations to try before giving up.
##         max.half(20 by default): the maximum number of times a step should be
##                                  halved before concluding that the step has
##                                  failed to improve the objective.
##         eps(1e-6 by default): the finite difference intervals to use when a
##                                Hessian function is not provided.
## Output: a list containing: f: the value of the objective function at the minimum.
##                            theta: the value of the parameters at the minimum.
##                            iter: the number of iterations taken to reach the minimum.
##                            g: the gradient vector at the minimum.
##                            Hi: the inverse of the Hessian matrix at the minimum.
  

  # approximate the hessian using the gradient if not provided
  if(is.null(hess)) {
    hess <- function(theta){
      grad_0 <- grad(theta, ...) # compute the gradient vector with the initial parameters
      n <- length(theta) # length of the parameter vector (n = # of parameters)
      hessian <- matrix(0, n, n) # finite difference Hessian (dimensions: n x n)
      for (i in 1:n){ # loop over parameters
        theta_1 <- theta; theta_1[i] <- theta_1[i] + eps # increase theta by eps
        grad_1 <- grad(theta_1, ...) # recompute the gradient vector with the perturbed parameters
        hessian[i,] <- (grad_1 - grad_0)/ eps # approximate second derivatives
      }
      hessian <- (t(hessian) + hessian) / 2 # making the hess matrix symmetric
      return(hessian)
    }
  }
  
  # Ensure all initial values are finite
  if (any(is.infinite(func(theta,...))) |
      any(is.infinite(grad(theta,...))) |
      any(is.infinite(hess(theta,...)))
      ){
    stop("The function or its derivates are not finite at the initial parameters.")
  }
  
  # Check if the optimization is within the tolerance given the initial guess
  if (all(abs(grad(theta,...)) < tol * (abs(func(theta, ...)) + fscale))){
    num_iter <- 0 # number of iterations taken to reach the minimum 
    
    # Check if the Hessian is positive definite at convergence. if not warn the user.
    if (min(eigen(hess(theta, ...))$values) <= 0){
      warning("The Hessian is not positive definite at convergence and therefore its inverse is not provided.")
      
      # Constructing output list
      out <- list('f'=func(theta, ...),'theta'=theta, 'iter'=num_iter,
                  'g'=grad(theta, ...), 'Hi'=NULL)
      
      return(out)
    } # exit the function
    
    # Constructing output list
    out <- list('f'=func(theta, ...),'theta'=theta, 'iter'=num_iter,
                'g'=grad(theta, ...), 'Hi'=chol2inv(chol(hess(theta, ...))))
    
    return(out)
  } # exit the function
  
  theta_k <- theta
  for (iter in 1:maxit) { # loop over the number of Newton iterations to try 
    
    # Checking that the Hessian is positive definite and perturbing it to be so f it isn't
    H <- hess(theta_k, ...)
    # flag if the hessian is positive definite
    posdef <- FALSE
    
    # Multiple of identity matrix used to perturb Hessian
    I <- diag(ncol(H)) # create an identity matrix same dim as matrix H
    i <- 0 # initialize counter for number of trials to perturb the hessian matrix
    multiple <- 0 # Initialize multiple (set = 0 to check if the initial hessian matrix is positive definite)
    
    # # Multiple of identity matrix scaled to the order of the Hessian, used to perturb Hessian
    # I <- diag(ncol(H)) * mean(diag(H))
    # k <- 0
    
    while(!posdef){
      
      if (i>10){
        # Give up if not successfully perturbed
        stop("Could not perturb H to be positive definite.")
      }
      
      # Try cholesky decomposition, will succeed if positive definite.
      error <- try(chol(H + multiple*I),silent=TRUE)
      
      # Checks to see if the cholesky decomposition failed
      if (inherits(error, "try-error")){ # H is not positive definite
        
        if (i==0){ # if the initial hessian matrix is not positive definite then
          # start perturbing by setting the multiplier to a small multiple
          # of Hessian matrix inf-norm
          multiple <- norm(H, c("I")) * 1e-6 
        }
        else{ # repeatedly multiply the multiplier by 10
          multiple <- multiple * 10
        }
      }
      
      else{ # H is positive definite
        posdef <- TRUE
        # If the initial hessian matrix is positive definite, (i=0 => multiple = 0), then it stays the same.
        # Otherwise we add a multiple of the identity that makes it positive definite.
        H <- H + multiple*I
      }
      
      i <- i + 1 # increase counter

      # if (k>10){
      #   # Give up if not successfully perturbed
      #   stop("Could not perturb H to be positive definite.")
      # }
      # # Try cholesky decomposition, will succeed if positive definite.
      # error <- try(chol(H + (10^k)*I),silent=TRUE)
      # # Checks to see if the cholesky decomposition failed
      # if (inherits(error, "try-error")){
      #   k <- k+1
      # }
      # else{
      #   posdef <- TRUE
      #   # If the initial hessian matrix is positive definite, k=0 and it stays the same.
      #   # Otherwise we add a multiple of the identity that makes it positive definite.
      #   H <- H + (10^k)*I
      # }
      
    }
    
    stepsize <- 1.0 # Initial step size
    for (step in 1:max.half){
      # estimate new parameters that would decrease the objective function
      
      # step = - inv(Hessian matrix) %*% gradient vector
      # - Hessian matrix * step = gradient vector
      # solve with Cholesky
      # compute Cholesky factor of Hessian matrix
      R <- chol(H) # avoid computing twice
      theta_k1 <- theta_k - stepsize * backsolve(R, forwardsolve(t(R), grad(theta_k, ...)))
      
      # if we went too far or the function is infinite, halve step size
      if (func(theta_k1,...) >= func(theta_k,...) | is.infinite(func(theta_k1, ...))){
        stepsize <- stepsize / 2
      }
    }
    
    if (func(theta_k1,...) >= func(theta_k,...)){
      warning("The step failed to reduce the objective function even after
              the maximum number of halvings.")
    }
    
    # Check if we are done, the optimization is within the tolerance
    if (all(abs(grad(theta_k1,...)) < tol * (abs(func(theta_k1, ...)) + fscale))){
      theta_k <- theta_k1
      num_iter <- iter # number of iterations taken to reach the minimum
      break
    }
    
    # move on to new step
    theta_k <- theta_k1
  }
  
  # At this point, the optimizer should have converged.
  # if we haven't yet converged, warn the user
  if (all(abs(grad(theta_k,...)) >= tol * (abs(func(theta_k, ...)) + fscale))){
    num_iter <- maxit
    warning("The optimizer has not converged despite reaching the maximum 
            number of iterations.")
  }
  
  # # Check if the Hessian is positive definite at convergence. if not warn the user
  # tryCatch(
  #   expr = {chol(hess(theta_k,...))} ,
  #   error = function(e){
  #     message("The Hessian is not positive definite at convergence.")
  #   }
  # )
  # 
  # Check if the Hessian is positive definite at convergence. if not warn the user.
  if (min(eigen(hess(theta_k, ...))$values) <= 0){
    warning("The Hessian is not positive definite at convergence and therefore its inverse is not provided.")
    
    # Constructing output list
    out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
                'g'=grad(theta_k, ...), 'Hi'= NULL)
    return(out)
  }

  # Constructing output list
  out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
              'g'=grad(theta_k, ...), 'Hi'=chol2inv(chol(hess(theta_k, ...))))
  
  return(out)
  
} ## newt

## End of code