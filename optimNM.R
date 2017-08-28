optim.NM <- function(objective.fn, params.init, iter.max=250, abs.tol=0.0001)
{
  # http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
  # params.init - initial/guess parameters
  # X.len - number of parameters, N
  # X.vert - a matrix of N+1 rows, and N columns, each row is a vertex of the simplex of N+1 vertices
  # objective.fn - the objective function to be minimized
  
  if (length(params.init) == 1)
  {
    params.init <- c(params.init, 0)
    objective.fn_1 <- objective.fn
    objective.fn <- objective.fn_1(x[1]) + 0*x[2] 
  }
  
  # define the search algorithm that will find the local optimmum given the initial simplex
  get.optim.NM <- function(X.vert, params.init, objective.fn, iter.max=250, abs.tol=0.0001)
  {
    # input dimension
    X.len <- length(params.init)
    
    # initialize controls before iterations of searching
    iter <- 0; not.converged <- 1; not.max.iter <- 1
    X.optim <- params.init; f_X.optim <- objective.fn(X.optim)
    
    # while loop, iterations
    while (not.converged & not.max.iter)
    {
      # get values at vertices
      f_X.vert <- apply(X = X.vert, MARGIN = 1, FUN = objective.fn) 
      
      # order ascending X.vert and f(X.vert), by f(X.vert)
      X.order <- sort(f_X.vert, index.return = TRUE)$ix
      X.vert <- X.vert[X.order, ]
      f_X.vert <- f_X.vert[X.order]
      
      # get centroid (mean on each dimension) of all points except the worst
      X.centr <- apply(X = X.vert[1:X.len, ], MARGIN = 2, FUN = mean)
      
      # get reflected point
      X.refl <- X.centr + coef.refl*(X.centr - X.vert[X.len+1, ])
      f_X.refl <- objective.fn(X.refl)
      
      if ((f_X.vert[1] <= f_X.refl) & (f_X.refl < f_X.vert[X.len]))
      { 
        # if the reflected point is better than the second worst, but not better than the best...
        # ... then obtain a new simplex by replacing the worst point with the reflected point
        X.vert[X.len+1, ] <- X.refl 
        
      } else if (f_X.refl < f_X.vert[1]) {
        
        # if the reflected point is the best point so far
        # ... then compute the expanded point
        X.expan <- X.centr + coef.expan*(X.centr - X.vert[X.len+1, ])
        f_X.expan <- objective.fn(X.expan)
        
        # ... if the expanded point is better than the reflected point
        if (f_X.expan < f_X.refl)
        {
          # ... then obtain a new simplex by replacing the worst point with the expanded point
          X.vert[X.len+1, ] <- X.expan
          
        } else {
          # ... else obtain a new simplex by replacing the worst point with the reflected point
          X.vert[X.len+1, ] <- X.refl
        }
        
      } else {
        
        # ... reflected point is not better than second worst
        # ... then compute the contracted point
        X.contr <- X.centr + coef.contr*(X.centr - X.vert[X.len+1, ])
        f_X.contr <- objective.fn(X.contr)
        
        # ... if the contracted point is better than the worst point
        if (f_X.contr < f_X.vert[X.len+1])
        {
          # ... then obtain a new simplex by replacing the worst point with the contracted point
          X.vert[X.len+1, ] <- X.contr
        } else {
          # ... shrink the simplex: X = X1 + coef.shrink(X-X1)
          X.vert <- sweep(coef.shrink*sweep(X.vert, 2, X.vert[1, ], FUN = "-"), 2, X.vert[1, ], FUN="+")
        }  
        
      }
      
      # get values at vertices
      f_X.vert <- apply(X = X.vert, MARGIN = 1, FUN = objective.fn) 
      
      # order asc X.vert and f(X.vert)
      X.order <- sort(f_X.vert, index.return = TRUE)$ix
      X.vert <- X.vert[X.order, ]
      f_X.vert <- f_X.vert[X.order]
      
      # update controls
      iter <- iter + 1 
      not.max.iter <- (iter < iter.max)*1
      not.converged <- (abs(objective.fn(X.vert[X.len, ])- objective.fn(X.vert[1, ])) > abs.tol)*1
      
      X.optim <- X.vert[1, ]; f_X.optim <- objective.fn(X.optim)
    }
    
    return(list(X.optim=X.optim, f_X.optim=f_X.optim))
    
  }
  
  # Nelder-Mead algorithm coefficients (as suggested in the original paper)
  coef.refl <- 1      # alpha
  coef.expan <- 2     # gamma
  coef.contr <- -0.5  # rho
  coef.shrink <- 0.5  # sigma
  
  # input dimension
  X.len <- length(params.init)
  
  # initialize a pocket of different starting vertices, so to minimize the risk of finding a local optimimum
  pocket.points <- list(lo = t(cbind(as.vector(params.init), params.init + diag(X.len)*(-0.50*params.init - 0.05))), 
                        zero = t(cbind(as.vector(params.init), params.init*0^diag(X.len))), 
                        up = t(cbind(as.vector(params.init), params.init + diag(X.len)*(0.50*params.init + 0.05)))) 
  
  # get optimums starting at each point in pocket
  optim.out <- get.optim.NM(X.vert=pocket.points$lo, params.init, objective.fn, iter.max, abs.tol)
  pocket <- c(optim.out$f_X.optim, optim.out$X.optim)
  optim.out <- get.optim.NM(X.vert=pocket.points$zero, params.init, objective.fn, iter.max, abs.tol)
  pocket <- rbind(pocket[], c(optim.out$f_X.optim, optim.out$X.optim))
  optim.out <- get.optim.NM(X.vert=pocket.points$up, params.init, objective.fn, iter.max, abs.tol)
  pocket <- rbind(pocket[], c(optim.out$f_X.optim, optim.out$X.optim))
  pocket <- pocket[order(pocket[, 1]), ]
  
  # re-initialize at best optimum in pocket, with a small simplex 
  optim.init <- pocket[1, 2:(X.len+1)]
  X.vert <- t(cbind(as.vector(optim.init), optim.init + diag(X.len)*(0.05*optim.init + 0.01)))
  optim.out <- get.optim.NM(X.vert=X.vert, params.init=optim.init, objective.fn, iter.max, abs.tol)
  
  X.optim <- list(par=optim.out$X.optim, value=optim.out$f_X.optim)
  
  return(X.optim)
}
