#' French swimmers performances data on 100m freestyle events
#'
#' A subset of data from reported performances of French swimmers during
#' 100m freestyle competitions between 2002 and 2016. See
#' https://link.springer.com/article/10.1007/s10994-022-06172-1 and
#' https://www.mdpi.com/2076-3417/8/10/1766 for dedicated description and
#' analysis.
#'
#' @format ## `swimmers`
#' A data frame with 76,832 rows and 4 columns:
#' \describe{
#'   \item{ID}{Indentifying number associated to each swimmer}
#'   \item{Input}{Age in years}
#'   \item{Output}{Performance in seconds on a 100m freestyle event}
#'   \item{Gender}{Competition gender}
#' }
#' @source <https://ffn.extranat.fr/webffn/competitions.php?idact=nat>
"swimmers"

#' Weight follow-up data of children in Singapore
#'
#' A subset of data from the GUSTO project (https://www.gusto.sg/) collecting
#' the weight over time of several children in Singapore.
#' See https://arxiv.org/abs/2011.07866 for dedicated description and
#' analysis.
#'
#' @format ## `weight`
#' A data frame with 3,629 rows and 4 columns:
#' \describe{
#'   \item{ID}{Indentifying number associated to each child}
#'   \item{sex}{Biological gender}
#'   \item{Input}{Age in months}
#'   \item{Output}{Weight in kilograms}
#' }
#' @source <https://www.gusto.sg/>
"weight"

