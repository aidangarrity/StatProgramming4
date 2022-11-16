

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,
                 fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
  # x_k+1 = x_k - [f''(x_k)]^-1 %*% f'(x_k)
  theta_k <- theta
  for (iter in 1:maxit) {
    
    theta_k1 <- theta_k - solve(hess(theta_k, ...)) %*% grad(theta_k, ...)
    
    if (all(abs(grad(theta_k,...)) < tol * (abs(func(theta_k, ...)) + fscale)) ){
      num_ter <- iter
      break
    }
    else{
      theta_k <- theta_k1
    }
  }

  out <- list('f'=func(theta_k, ...),'theta'=theta_k, 'iter'=num_iter,
              'g'=grad(theta_k, ...), 'Hi'=solve(hess(theta_k, ...)))
  
  
}
