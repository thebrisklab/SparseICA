# Sparse independent component analysis (ICA)

## 1. Data

You can directly use `X.csv` and `S.csv` as the generated observed data and true signal data. Or you can do it by yourself after loading the `jngcaFunctions.R`. Thanks to Benjamin Risk （<brisk@emory.edu>） for the functions.  

## 2. Refined ICA

The `refinedICA.R` contains the main part of our refined ICA algorithm. There are four different algorithms allowed in one function: `refinedICA()`. 

You can have the following choices:  

+ Relax logistic: this one executed a little bit slowly, especially when compared with others. 
+ Relax Softmax: it can run very fast, and it allows sparsity results. 
+ FastICA logistic We directly embedded the codes of `FastICA`  for it. So it is just to call the `FastICA` codes with logistic density. 
+ FastICA tanh: This is also the `FastICA`  codes, under tanh function. 

Note: For *relax_laplace* and *relax_logistic* algorithms, you can check the convergence situation by using `converge_plot = T`. However, there is no convergence plot for FastICA methods, since the embedded codes didn't allow the convergence plot. 

