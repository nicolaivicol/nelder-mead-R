# Nelder-Mead simplex optimization method in R
The code in R for Nelder–Mead simplex method to find the minimum of an objective function in a multidimensional space    
Read more about Nelder-Mead method: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method    

# Instructions
The function `optim.NM` returns the minimum of a multidimensional objective function:   
```R
optim.NM(objective.fn, params.init, iter.max, abs.tol)
```
with the standard inputs for European and American options:       
- `objective.fn`:  objective function to be minimized
- `params.init`:  initial guess for parameters of the function, a vector of zero's can be provided
- `iter.max`:  maximum iterations fo the search
- `abs.tol`:  tolerated precision of the objective function value when the search can be stopped

# Examples
The Rosenbrock function with `a=1`, `b=100` has the global minimum at `X=(a, a^2)`.   
https://en.wikipedia.org/wiki/Rosenbrock_function
```R
obj.fn.rosenbrock <- function(X)
{
  a <- 1
  b <- 100
  d <- length(X)
  Xi <- X[1:(d-1)]
  Xnext <- X[2:d]
  y <- sum(b*(Xnext-Xi^2)^2 + (Xi-a)^2)
  return(y)
}
```
```R
optim.NM(udf.rosenbrock, c(0, 0), 1000, 10^-9)
```
```R
$par
[1] 0.9999877 0.9999738

$value
[1] 3.957329e-10
```
```R
optim.NM(obj.fn.rosenbrock, c(0, 0, 0), 1000, 10^-9)
```
```R
$par
[1] 1.000004 1.000007 1.000015

$value
[1] 3.165469e-10
```
