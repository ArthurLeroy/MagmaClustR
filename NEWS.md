# MagmaClustR (development version)

## Major
* Remove the package 'optimr' dependency and switch to base 'optim()' function
* Increase convergence tolerance in 'optim()', which was too slow

## Minor
* Remove error message in 'train_magmaclust()' when common_hp_k = FALSE
* Change the default intervals for hyper-parameters in 'simu_db()'
* Automatically remove rows with missing data
* Change position of the 'grid_inputs' argument in prediction functions
* Remove the internal functions from the index documentation 
* Fix 'ID' in hyperposterior() and hyperposterior_clust() when not character

# MagmaClustR 1.0.0
* Initial release
