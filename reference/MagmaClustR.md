# MagmaClustR : Clustering and Prediction using Multi-Task Gaussian Processes

The **MagmaClustR** package implements two main algorithms, called
*Magma* and *MagmaClust*, using a multi-task GPs model to perform
predictions for supervised learning problems. Theses approaches leverage
the learning of cluster-specific mean processes, which are common across
similar tasks, to provide enhanced prediction performances (even far
from data) at a linear computational cost (in the number of tasks).
*MagmaClust* is a generalisation of *Magma* where the tasks are
simultaneously clustered into groups, each being associated to a
specific mean process. User-oriented functions in the package are
decomposed into training, prediction and plotting functions. Some basic
features of standard GPs are also implemented.

## Details

For a quick introduction to MagmaClustR, please refer to the README at
<https://github.com/ArthurLeroy/MagmaClustR>

## Author(s)

Arthur Leroy, Pierre Pathe and Pierre Latouche  
Maintainer: Arthur Leroy - <arthur.leroy.pro@gmail.com>

## References

Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey.  
MAGMA: Inference and Prediction with Multi-Task Gaussian Processes.
*Machine Learning*, 2022,
<https://link.springer.com/article/10.1007/s10994-022-06172-1>

Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey.  
Cluster-Specific Predictions with Multi-Task Gaussian Processes.
*Journal of Machine Learning Research*, 2023,
<https://jmlr.org/papers/v24/20-1321.html>

## Examples

### Simulate a dataset, train and predict with Magma 

set.seed(4242)  
data_magma \<- simu_db(M = 11, N = 10, K = 1)  
magma_train \<- data_magma %\>% subset(ID %in% 1:10)  
magma_test \<- data_magma %\>% subset(ID == 11) %\>% head(7)  

magma_model \<- train_magma(data = magma_train)  
magma_pred \<- pred_magma(data = magma_test, trained_model =
magma_model, grid_inputs = seq(0, 10, 0.01))  

### Simulate a dataset, train and predict with MagmaClust 

set.seed(4242)  
data_magmaclust \<- simu_db(M = 4, N = 10, K = 3)  
list_ID = unique(data_magmaclust\$ID)  
magmaclust_train \<- data_magmaclust %\>% subset(ID %in%
list_ID\[1:11\])  
magmaclust_test \<- data_magmaclust %\>% subset(ID == list_ID\[12\])
%\>% head(5)  

magmaclust_model \<- train_magmaclust(data = magmaclust_train)  
magmaclust_pred \<- pred_magmaclust(data = magmaclust_test,  
trained_model = magmaclust_model, grid_inputs = seq(0, 10, 0.01))  

## See also

Useful links:

- <https://github.com/ArthurLeroy/MagmaClustR>

- <https://arthurleroy.github.io/MagmaClustR/>

- Report bugs at <https://github.com/ArthurLeroy/MagmaClustR/issues>

## Author

**Maintainer**: Arthur Leroy <arthur.leroy.pro@gmail.com>
([ORCID](https://orcid.org/0000-0003-0806-8934))

Authors:

- Pierre Latouche <pierre.latouche@gmail.com>

Other contributors:

- Pierre Pathé <pathepierre@gmail.com> \[contributor\]

- Alexia Grenouillat <grenouil@insa-toulouse.fr> \[contributor\]

- Hugo Lelievre <lelievre@insa-toulouse.fr> \[contributor\]
