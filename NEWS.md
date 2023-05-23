# MagmaClustR (development version)

# MagmaClustR 1.1.2

## Minor 
* Fix a bug occurring in pred_magmaclust() for a 'trained_model' with hp_i = FALSE
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
