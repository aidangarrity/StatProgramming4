## Aidan Garrity (s1997567), Eleni Michaelidou (s2226022)

## Address of Githup repo: https://github.com/aidangarrity/StatProgramming4.git

## Brief description of what each team member contributed to the project:

## Proportion of the work was undertaken by each team member: 
## Aidan Garrity: 
## Eleni Michaelidou:


## ------------------------- Overview of the code ------------------------------


## -----------------------------------------------------------------------------

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,
                 fscale=1,maxit=100,max.half=20,eps=1e-6) {
## Function for ...
## Inputs: theta: a vector of initial values for the optimization parameters.
##         func: the objective function to minimize.
##         grad: the gradient function.
##         hess(NULL by default): the Hessian matrix function (if not supplied
##                                then an approximation to the Hessian by finite
##                                differencing of the gradient vector is obtained).
##         toll(1e-8 by default): the the convergence tolerance.
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
  
<<<<<<< HEAD
  ## approximate the hessian using the gradient if not provided
  if(is.null(hess)) {
    hess <- function(theta){
      grad_0 <- grad(theta, ...) ## compute the gradient vector with the initial parameters
      n <- length(theta) ## length of the parameter vector (n = # of parameters)
      hessian <- matrix(0, n, n) ## finite difference Hessian (dimensions: n x n)
      for (i in 1:n){ ## loop over parameters
        theta_1 <- theta; theta_1[i] <- theta_1[i] + eps ## increase theta by eps
        grad_1 <- grad(theta_1, ...) ## recompute the gradient vector with the perturbed parameters
        hessian[i,] <- (grad_1 - grad_0)/ eps ## approximate second derivatives
      }
      hessian <- (t(hessian) + hessian) / 2 ## making the hess matrix symmetric
      return(hessian)
    }
    
=======
  # approximate the hessian using the gradient if not provided
  if(is.null(hess)){
    grad_0 <- grad(theta, ...) # compute the gradient vector with the initial parameters
    n <- length(theta) # length of the parameter vector (n = # of parameters)
    hess <- matrix(0, n, n) # finite difference Hessian (dimensions: n x n)
    for (i in 1:n){ # loop over parameters
      theta_1 <- theta; theta_1[i] <- theta_1[i] + eps # increase theta by eps
      grad_1 <- grad(theta_1, ...) # recompute the gradient vector with the perturbed parameters
      hess[i,] <- (grad_1 - grad_0)/ eps # approximate second derivatives
    }
    
    hess <- (t(hess) + hess) / 2 # making the hess matrix symmetric
>>>>>>> 101e02c26466d138ead7e92124d349ca294b54b8
  }
  
  # Ensure all initial values are finite
  if (any(is.infinite(func(theta,...))) |
      any(is.infinite(grad(theta,...))) |
      any(is.infinite(hess(theta,...)))
      ){
    stop("The function or its derivates are not finite at the initial parameters.")
  }
  
  
  # x_k+1 = x_k - [f''(x_k)]^-1 %*% f'(x_k)
  theta_k <- theta
  for (iter in 1:maxit) { # loop over the number of Newton iterations to try 
    
    stepsize <- 1.0
    for (step in 1:max.half){
      # estimate new parameters that would decrease the objective function
<<<<<<< HEAD
      theta_k1 <- theta_k - stepsize* (chol2inv(chol(hess(theta_k, ...))) %*% grad(theta_k, ...))
=======
      
      # step = - inv(Hessian matrix) %*% gradient vector
      # - Hessian matrix * step = gradient vector
      
      # solve with Cholesky
      # compute Cholesky factor of Hessian matrix
      R <- chol(hess(theta_k, ...)) # avoid computing twice
      theta_k1 <- theta_k - stepsize * backsolve(R, forwardsolve(t(R), grad(theta_k, ...)))
      
      #theta_k1 <- theta_k - stepsize* (solve(hess(theta_k, ...)) %*% grad(theta_k, ...))
>>>>>>> 101e02c26466d138ead7e92124d349ca294b54b8
      
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
  
<<<<<<< HEAD
  tryCatch(
    expr = {chol(hess(theta_k,...))} ,
    error = function(e){
      message("The Hessian is not positive definite at convergence.")
    }
  )
=======
  try(chol(hess(theta_k, ...)))
  
  # Check if the Hessian is positive definite at convergence.
  # if not warn the user.
  if (min(eigen(hess(theta_k, ...))$values) <= 0){
    warning("The Hessian is not positive definite at convergence.")
  }
>>>>>>> 101e02c26466d138ead7e92124d349ca294b54b8

  out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
              'g'=grad(theta_k, ...), 'Hi'=chol2inv(chol(hess(theta_k, ...))))
  
  
} ## newt

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}
print(newt(c(2,-2), rb, gb, hb, 2))


ra <- function(th){
  2*th[1]^2 + th[2]^2-4 + 6*th[2] - 7*th[1]
}
ga <- function(th){
  c(4*th[1]-7, 2*th[2]+6)
}
ha <- function(th){
  h <- matrix(0,2,2)
  h[1,1] <- 4
  h[2,2] <- 2
  h[1,2] <- 0
  h[2,1] <- 0
  h
}

print(newt(c(100,100), ra,ga,ha))
