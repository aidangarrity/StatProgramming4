

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,
                 fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
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
    
    # Check if we are done, the opimization is within the tolerance
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
  
  
}

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
