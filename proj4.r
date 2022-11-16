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
  
  # approximate the hessian using the gradient if not provided
  if(is.null(hess)){
    
  }
  
  # Ensure all initial values are finite
  if (any(is.infinite(func(theta,...))) |
      any(is.infinite(grad(theta,...))) |
      any(is.infinite(hess(theta,...)))){
    stop("The function or its derivates are not finite at the initial parameters.")
  }
  
  
  # x_k+1 = x_k - [f''(x_k)]^-1 %*% f'(x_k)
  theta_k <- theta
  for (iter in 1:maxit) {
    
    stepsize <- 1.0
    for (step in 1:max.half){
      # estimate new parameters that would decrease the objective function
      theta_k1 <- theta_k - stepsize* (solve(hess(theta_k, ...)) %*% grad(theta_k, ...))
      
      # if we went to far, halve step size
      if (func(theta_k1,...) >= func(theta_k,...)){
        stepsize <- stepsize / 2
      }
    }
    
    if (func(theta_k1,...) >= func(theta_k,...)){
      warning("The step failed to reduce the objective function even after
              the maximum number of halvings.")
    }
    
    # Check if we are done, the optimization is within the tolerance
    if (all(abs(grad(theta_k,...)) < tol * (abs(func(theta_k, ...)) + fscale)) ){
      num_ter <- iter
      break
    }
    
    # move on to new step
    else{
      theta_k <- theta_k1
    }
    
  }
  
  # At this point, the opimizer should have converged.
  # if we haven't yet converged, warn the user
  if (all(abs(grad(theta_k,...)) >= tol * (abs(func(theta_k, ...)) + fscale))){
    warning("The optimzer has not converged despite reaching the maximum 
            number of iterations.")
  }
  
  try(chol(hess(theta_k, ...)))

  out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
              'g'=grad(theta_k, ...), 'Hi'=solve(hess(theta_k, ...)))
  
  
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
