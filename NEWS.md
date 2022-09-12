# MagmaClustR (development version)

## Major
* Remove the package 'optimr' dependency and switch to base 'optim()' function
* Increase convergence tolerance in 'optim()', which was too slow
* Implement expand_grid_inputs() to create an n-dimensional
grid on which to evaluate the GP.
* Implement regularize_data() to approximate data on a grid,
and control the size of the covariance matrix.
* Add a Reference column to control how the algorithm works with
multidimensional data.
* Change simu_db() to have more interesting 2D data.

## Minor
* Remove error message in 'train_magmaclust()' when common_hp_k = FALSE
* Change the default intervals for hyper-parameters in 'simu_db()'
* Automatically remove rows with missing data
* Change position of the 'grid_inputs' argument in prediction functions
* Remove the internal functions from the index documentation 
* Fix 'ID' in hyperposterior() and hyperposterior_clust() when not character
* Round inputs to 6 significant digits to avoid numerical errors.
* Generalize the creation of a grid in any dimension in case grid_inputs is not
specified in the prediction functions.
* Use of match.closest() from the MALDIquant package.

# MagmaClustR 1.1.0
* First major update
