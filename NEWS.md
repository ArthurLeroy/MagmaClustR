# MagmaClustR (development version)

## Major

## Minor

* Fix issues when plotting in 2D mean process predictions 

# MagmaClustR 1.2.1

## Major
* Make 'get_full_cov = T' the default option when using pred_magmaclust()

## Minor
* Remove the useless messages when using samples=T in plotting functions
* Create an error when users try train_magmaclust(nb_cluster = 1)
* Make 'samples = T' the default plot when using prediction functions
* Add mention of the 'hp()' function to the documentation
* Fix Warning message when plotting clusters with exact same posterior probability
* Display plots for the generic (no data provided) prediction that are condistent with classic ones
* Disable the symmetry check when generating posterior samples for predictions

# MagmaClustR 1.2.0

## Major

* Add new functions sample_gp() and sample_magmaclust() to sample from the posterior predictions of GPs, Magma, and MagmaClust
* Propose a new visualisation based on posterior samples instead of Credible Intervals for both plot_gp() and plot_magmaclust()
* Allow pred_*() functions to be used without providing the 'data' argument by returning the mean processes (or priors) as generic predictions. 

## Minor 

* Add an option to generate multiple curves from sample_gp()
* Fix a bug when the 'Reference' column is present when using plot_gp() 
* Fix the '\docType{package}' bug in roxygen2
* Add 'plot_magma()' as a duplicated name for 'plot_gp()'
* Remove the useless 'Reference' column in the prediction of the mean processes
* Fix a bug when both 'pred' and 'samples' are provided in plot_samples()

# MagmaClustR 1.1.2

## Minor 
* Fix a bug occurring in pred_magmaclust() for a 'trained_model' with hp_i=FALSE
* Simplify the use of hyperposterior() and hyperposterior_clust() by providing a 'trained_model' argument.

# MagmaClustR 1.1.1

## Minor
* Fix an issue regarding deprecation of .data in 'tidyselect'
* Fix a unit testing issue

# MagmaClustR 1.1.0

## Major
* Provide 4 vignettes explaining in detail how the different features of MagmaClustR work on practical examples. 
* Implement expand_grid_inputs() to help create customised n-dimensional input
grids on which to evaluate the GP.
* Implement regularize_data() to project a dataset on a specific input grid,
(possibly to control the size of the resulting covariance matrices and the associated running time).
* Add an internal 'Reference' column to datasets, to provide an adequate identifier for multidimensional inputs.
* Implement a new version of simu_db() to generate more realistic 2-D datasets.

## Minor
* Round inputs to 6 significant digits to avoid numerical errors.
* Generalise the creation of a grid in any dimension when 'grid_inputs' is not
specified in the prediction functions.

# MagmaClustR 1.0.1

## Major
*Remove the package 'optimr' dependency and switch to base 'optim()' function
*Increase convergence tolerance in 'optim()', which was too slow

## Minor
*Fix the warnings about the absolute value function in the Cpp code
*Remove error message in 'train_magmaclust()' when common_hp_k = FALSE
*Change the default intervals for hyper-parameters in 'simu_db()'
*Automatically remove rows with missing data
*Change position of the 'grid_inputs' argument in prediction functions
*Remove the internal functions from the index documentation
*Fix 'ID' in hyperposterior() and hyperposterior_clust() when not character


# MagmaClustR 1.0.0
Initial release
