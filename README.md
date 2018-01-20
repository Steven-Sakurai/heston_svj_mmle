# Heston SVJ MLE  

Written with CMake.    
In progress.  

### Task List:  
1. Calculate the initial truncated filter moments. Done  
2. Lagrange cubic interpolation coef. Done  
3. Make fixes due to jump. Done  
4. Define Functors to allow for explicit calculations. Done  
5. Calculate `MarTranCoef` `alpha` Done(not consistent with previous matlab code, need confirmation)    
6. Calculate `BasisCoef` matrix.  
7. Implement induction algorithm.     
8. Obtain the function l(\theta).  
9. Maximize the log-likelihood function, find the estimators.    