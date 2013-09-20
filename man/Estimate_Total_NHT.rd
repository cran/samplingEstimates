\name{Estimate.Total.NHT}
\alias{Estimate.Total.NHT}
\title{Narain-Horvitz-Thompson estimates for a total from survey data}
\description{Produces estimates for population totals that are estimated using the Narain (1951); Horvitz-Thompson (1952) point estimator and from survey data obtained from a single-stage sampling design, i.e. direct element sampling.  }
\usage{
Estimate.Total.NHT(MatY.s                ,
                   VecWk.s               ,
                   VarEst         = "SYG",
                   MatPkl.s       = NULL ,
                   PopSize        = NULL ,
                   VecStratLb.s   = NULL ,
                   VecStratSize.s = NULL ,
                   ShowStrata     = FALSE,
                   VecDomainLb.s  = NULL )
}
\arguments{
\item{MatY.s}{matrix (dataframe or vector) with \eqn{n} rows (observations) and \eqn{Q} columns (variables of interest), where \eqn{n} is the overall sample size of elements. The argument \code{MatY.s} can also be a vector which is internally treated as a matrix with \eqn{Q=1}. There must not be any missing value.}
\item{VecWk.s}{vector of the elements sampling weights; its length is equal to \eqn{n}, the sample size. Values in \code{VecWk.s} must be greater than or equal to one. Columns of \code{MatY.s} and length of \code{VecWk.s} must be the same. There must not be any missing value.}
\item{VarEst}{string indicating the mathematical expression for estimating the variance. Available options are: \code{"HT"}, \code{"SYG"} or \code{"Hajek"}. If \code{VarEst} argument is omitted, the default is \code{"SYG"}, which requires the provision of the matrix of joint inclusion probabilities \code{MatPkl.s}. In the case that the argument \code{MatPkl.s} is not provided, \code{VarEst} is set to \code{"Hajek"}. When using \code{VarEst="Hajek"}, it is assumed a high-entropy sampling design, see Hajek (1964); care should be taken with highly-stratified samples, e.g. Berger (2005).}
\item{MatPkl.s}{matrix of the second-order inclusion probabilities; its number of rows and columns is equal to \eqn{n}, the overall sample size of observed elements. Values in \code{MatPkl.s} must be greater than zero and less than or equal to one. There must not be any missing value.}
\item{PopSize}{population size \eqn{N}. This argument may be optional; if it is not provided the computations are made using \eqn{\hat{N}=\sum_{k\in s}w_k}, which estimates the total of elements in the population.}
\item{VecStratLb.s}{vector of the strata labels; its length is equal to \eqn{n}, the sample size. Values in the argument \code{VecStratLb.s} can be numeric (integers), strings or a factor. It does not need to be sorted, however, all other arguments (variables, vectors, matrices) of size \eqn{n} must follow the same order of \code{VecStratLb.s} correspondingly. This argument is optional, if it is not provided the computations are made assuming that there is no stratification. There must not be any missing value.}
\item{VecStratSize.s}{vector of the strata population sizes; its length is equal to \eqn{n}, the sample size. This vector contains, for each of the \eqn{n} observations, the size of the stratum each observation belongs to. The argument \code{VecStratSize.s} does not need to be sorted, however, all other arguments (variables, vectors, matrices) of size \eqn{n} must follow the same order of \code{VecStratSize.s} correspondingly. There must not be any missing value.}
\item{ShowStrata}{logical. If \code{TRUE} partial results from each stratum is displayed. This is an optional argument; default is \code{FALSE}.}
\item{VecDomainLb.s}{vector of the domains (sub-groups) labels; its length is equal to \eqn{n}, the sample size. Values in the argument \code{VecDomainLb.s} can be numeric (integers) or strings (characters). They do not need to be sorted, however this variable (column) must follow the same order of \code{VecStratLb.s} correspondingly with the other variables, vectors or matrices of size \eqn{n}. This argument is optional and there must not be any missing value.}
}
\details{
For the population total of the variable \eqn{y}:
\deqn{t = \sum_{k\in U} y_k}
the unbiased Narain (1951); Horvitz-Thompson (1952) estimator of \eqn{t} is given by:
\deqn{\hat{t}_{NHT} = \sum_{k\in s} w_k y_k}
where \eqn{w_k} denotes the sampling weight of the \eqn{k}-th element in the sample \eqn{s}, \eqn{w_k=1/\pi_k} with \eqn{\pi_k} denoting the inclusion probability of the \eqn{k}-th element in the sample. Let \eqn{\pi_{kl}} denotes the joint-inclusion probabilities of the \eqn{k}-th and \eqn{l}-th elements in the sample \eqn{s}. The variance of \eqn{\hat{t}_{NHT}} is given by:
\deqn{V(\hat{t}_{HT}) = \sum_{k\in U}\sum_{l\in U} (\pi_{kl}-\pi_k\pi_l)w_k y_k w_l y_l}
which can therefore be estimated by the Horvitz-Thompson variance estimator (implemented by the current function if \code{VarEst="HT"}):
\deqn{\hat{V}_{HT}(\hat{t}_{NHT}) = \sum_{k\in s}\sum_{l\in s} \frac{\pi_{kl}-\pi_k\pi_l}{\pi_{kl}}w_k y_k w_l y_l}

If the utilised sampling design is of fixed-size, the variance \eqn{V(\hat{t}_{NHT})} can be estimated by the Sen-Yates-Grundy variance estimator (implemented by the current function if \code{VarEst="SYG"}):
\deqn{\hat{V}_{SYG}(\hat{t}_{NHT}) = \frac{-1}{2}\sum_{k\in s}\sum_{l\in s} \frac{\pi_{kl}-\pi_k\pi_l}{\pi_{kl}}\left(w_k y_k - w_l y_l\right)^2}
For large-entropy sampling designs, the variance of \eqn{\hat{t}_{NHT}} is approximated by the Hajek (1964) variance:
\deqn{V_{Hajek}(\hat{t}_{NHT}) \doteq \frac{N}{N-1}\left[\sum_{k\in U} w_k y_k^2\left(\frac{w_k-1}{w_k}\right)-dG^2\right]}
with \eqn{d=\sum_{k\in U}w_k^{-2}(w_k-1)} and \eqn{G=d^{-1}\sum_{k\in U}w_k^{-1}(w_k-1)y_k}.

This approximate variance can therefore be estimated by the variance estimator (implemented by the current function if \code{VarEst="Hajek"}):
\deqn{\hat{V}_{Hajek}(\hat{t}_{NHT}) = \frac{n}{n-1}\left[\sum_{k\in s}w_k^2 y_k^2\left(\frac{w_k-1}{w_k}\right)-\hat{d}\hat{G}^2\right]}
where \eqn{\hat{d}=\sum_{k\in s}w_k^{-1}(w_k-1)} and \eqn{\hat{G}=\hat{d}^{-1}\sum_{k\in s}(w_k-1)y_k}.

The Hajek (1964) variance approximation is designed for large-entropy sampling designs and large populations, i.e. care should be taken with highly-stratified samples, e.g. Berger (2005).
  }
\value{
The function returns a dataframe with \eqn{Q} rows (the number of variables of interest) and some columns depending on input information and used expressions in computations. The results in the returned columns are:
\item{Statistic}{the utilised point estimator.}
\item{VariableName}{the name of the variable of interest.}
\item{Estimate}{the point estimate obtained from evaluating the sample data.}
\item{Variance}{the estimated variance of the point estimator.}
\item{StdErr}{the estimated standard error of the point estimator.}
\item{AbsErr}{the estimated absolute error of the point estimator.}
\item{LInfCI95}{the lower limit of the 95 percent confidence interval.}
\item{LSupCI95}{the upper limit of the 95 percent confidence interval.}
\item{Range95}{the range (width) of the 95 percent confidence interval.}
\item{PctCVE}{the estimated coefficient of variation (in percentage).}
\item{DEff}{the estimated design effect.}
\item{n}{the overall sample size.}
\item{Nhat}{an estimate of the population size (number of elements in the population) \eqn{\hat{N}=\sum_{k\in s}w_k}.}
\item{fhat}{an estimate of the overall sampling fraction \eqn{\hat{f}=n/\hat{N}}.}
\item{N}{the population size (total of elements in the population).}
\item{f}{the overall sampling fraction.}
If a stratified sampling design was specified and if \code{ShowStrata=TRUE} some further columns are displayed with partial results. Note that these per-stratum partial results are NOT returned by the function, they are only on-screen information.
\item{h}{stratum counter.}
\item{Stratum}{stratum label (integer, character).}
\item{nh}{the sample size for the stratum \eqn{h}.}
\item{Nh}{the size of the stratum \eqn{h} (total of elements in the stratum \eqn{h}).}
\item{fh}{the sampling fraction for the stratum \eqn{h}.}
\item{Wh}{the relative weight of the stratum \eqn{h} among all strata \eqn{W_{h}=n_{h}/N_{h}}.}
If domains of study were specified these extra columns are displayed. Note that these per-domain results are NOT returned by the function, they are only on-screen information.
\item{d}{domain counter.}
\item{Domain}{domain label.}
\item{nd}{the sample size in the domain \eqn{d}.}
\item{Ndhat}{an estimate of the population size (number of elements) for the domain \eqn{d}.}
\item{fdhat}{an estimate of the sampling fraction for the domain \eqn{d}.}
\item{Wdhat}{an estimate of the relative weight of the domain \eqn{d} among all domains.}
}
\references{
Berger, Y. G. (2005) Variance estimation with highly stratified sampling designs with unequal probabilities. \emph{Australian & New Zealand Journal of Statistics}, \bold{47}, 365--373.

Hajek, J. (1964) Asymptotic theory of rejective sampling with varying probabilities from a finite population. \emph{The Annals of Mathematical Statistics}, \bold{35}, 4, 1491--1523.

Horvitz, D. G. and Thompson, D. J. (1952) A generalization of sampling without replacement from a finite universe. \emph{Journal of the American Statistical Association}, \bold{47}, 663--685.

Narain, R. D. (1951) On sampling without replacement with varying probabilities. \emph{Journal of the Indian Society of Agricultural Statistics}, \bold{3}, 169--175.

Sen, A. R. (1953) On the estimate of the variance in sampling with varying probabilities. \emph{Journal of the Indian Society of Agricultural Statistics}, \bold{5}, 119--127.

Yates, F. and Grundy, P. M. (1953) Selection without replacement from within strata with probability proportional to size. \emph{Journal of the Royal Statistical Society B}, \bold{15}, 253--261.
}

\examples{
##################################
## Setting up data to run examples
##################################
data(Sample1)           ## Loads a data frame with the sample to be used in examples
N         <- 570        ## Defining the population size
## Approximating the 2nd order inclusion probabilities with sample based quantitites
## (Note: this approximation is only suitable for large-entropy sampling designs)
require(samplingVarEst) ## Loading the necessary package
Probs2Mat <- Pkl.Hajek.s(Sample1$InclProbs) ## function from samplingVarEst package
head(Sample1)           ## Showing the first rows of the sample data to be used

############################################################
## Example 1: A variable of interest, without stratification
############################################################
Estimate.Total.NHT(MatY.s   = Sample1$y1     ,
                   VecWk.s  = Sample1$Weights)

Estimate.Total.NHT(MatY.s   = Sample1$y1     ,
                   VecWk.s  = Sample1$Weights,
                   VarEst   = "HT"           )

Estimate.Total.NHT(MatY.s   = Sample1$y1     ,
                   VecWk.s  = Sample1$Weights,
                   VarEst   = "SYG"          ,
                   MatPkl.s = Probs2Mat      )

Estimate.Total.NHT(MatY.s   = Sample1$y1     ,
                   VecWk.s  = Sample1$Weights,
                   VarEst   = "SYG"          ,
                   MatPkl.s = Probs2Mat      ,
                   PopSize  = N              )


###################################################################################
## Example 2: A matrix/dataframe of 2 variables of interest, without stratification
###################################################################################
Estimate.Total.NHT(MatY.s   = Sample1[ ,c("y1","y2")],
                   VecWk.s  = Sample1$Weights        ,
                   VarEst   = "SYG"                  ,
                   MatPkl.s = Probs2Mat              ,
                   PopSize  = N                      )


#########################################################
## Example 3: A variable of interest, with stratification
#########################################################
Estimate.Total.NHT(MatY.s         = Sample1$y1             ,
                   VecWk.s        = Sample1$Weights        ,
                   VecStratLb.s   = Sample1$CharStrataNames,
                   VecStratSize.s = Sample1$StrataSizes    )

Estimate.Total.NHT(MatY.s         = Sample1$y1             ,
                   VecWk.s        = Sample1$Weights        ,
                   VecStratLb.s   = Sample1$CharStrataNames,
                   VecStratSize.s = Sample1$StrataSizes    ,
                   ShowStrata     = TRUE                   )


###############################################################################
## Example 4: A matrix/dataframe (2 variables of interest), with stratification
###############################################################################
Estimate.Total.NHT(MatY.s         = Sample1[ ,c("y1","y2")],
                   VecWk.s        = Sample1$Weights        ,
                   VecStratLb.s   = Sample1$CharStrataNames,
                   VecStratSize.s = Sample1$StrataSizes    ,
                   ShowStrata     = TRUE                   )


#################################################################################
## Example 5: A matrix/dataframe (2 variables), no strata, with unplanned domains
#################################################################################
Estimate.Total.NHT(MatY.s        = Sample1[ ,c("y1","y2")],
                   VecWk.s       = Sample1$Weights        ,
                   VecDomainLb.s = Sample1$CharDoms       )

Estimate.Total.NHT(MatY.s        = Sample1[ ,c("y1","y2")],
                   VecWk.s       = Sample1$Weights        ,
                   VecDomainLb.s = Sample1$NumDoms        )


###################################################################################
## Example 6: A matrix/dataframe (2 variables), with strata, with unplanned domains
###################################################################################
Estimate.Total.NHT(MatY.s         = Sample1[ ,c("y1","y2")],
                   VecWk.s        = Sample1$Weights        ,
                   VecStratLb.s   = Sample1$CharStrataNames,
                   VecStratSize.s = Sample1$StrataSizes    ,
                   ShowStrata     = TRUE                   ,
                   VecDomainLb.s  = Sample1$CharDoms       )
}
\keyword{estimates}
\keyword{total}
\keyword{strata}
\keyword{domains}
