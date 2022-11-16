

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
    for (step in max.half){
      theta_k1 <- theta_k - stepsize* (solve(hess(theta_k, ...)) %*% grad(theta_k, ...))
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

  out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
              'g'=grad(theta_k, ...), 'Hi'=solve(hess(theta_k, ...)))
  
  
}
