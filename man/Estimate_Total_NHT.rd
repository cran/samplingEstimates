\name{Estimate.Total.NHT}
\alias{Estimate.Total.NHT}
\title{Narain-Horvitz-Thompson estimates for a total from survey data}
\description{Produces estimates for population totals that are estimated using the Narain (1951); Horvitz-Thompson (1952) point estimator and from survey data obtained from a single stage sampling design, i.e. direct element sampling.  }
\usage{
Estimate.Total.NHT(MatY.s                  ,
                   VecWk.s                 ,
                   VarEst           = "SYG",
                   MatPkl.s         = NULL ,
                   PopSize          = NULL ,
                   VecStrataId.s    = NULL ,
                   VecStrataSizes.H = NULL ,
                   ShowStrata       = FALSE)
}
\arguments{
\item{MatY.s}{matrix (dataframe or vector) with \eqn{n} rows (observations) and \eqn{Q} columns (variables of interest), where \eqn{n} is the overall sample size of elements. The argument \code{MatY.s} can also be a vector which is internally treated as a matrix with \eqn{Q=1}. There must not be any missing value.}
\item{VecWk.s}{vector of the elements sampling weights; its length is equal to \eqn{n}, the sample size. Values in \code{VecWk.s} must be greater than or equal to one. Columns of \code{MatY.s} and length of \code{VecWk.s} must be the same. There must not be any missing value.}
\item{VarEst}{string indicating the mathematical expression for estimating the variance. Available options are: \code{"HT"}, \code{"SYG"} or \code{"Hajek"}. If \code{VarEst} argument is omitted, the default is \code{"SYG"}, which requires the provision of the matrix of joint inclusion probabilities \code{MatPkl.s}. In the case that the argument \code{MatPkl.s} is not provided, \code{VarEst} is set to \code{"Hajek"}. Here it is pertinent to note that When using \code{VarEst="Hajek"}, it is assumed that a high-entropy sampling design was utilised when drawing the sample, see Hajek(1964); care should be taken with highly-stratified samples, e.g. Berger (2005).}
\item{MatPkl.s}{matrix of the second-order inclusion probabilities; its number of rows and columns is equal to \eqn{n}, the overall sample size of observed elements. Values in \code{MatPkl.s} must be greater than zero and less than or equal to one. There must not be any missing value.}
\item{PopSize}{population size \eqn{N}. This argument may be optional; if it is not provided the computations are made using \eqn{\hat{N}=\sum_{k\in s}w_k}, which is itself the Narain-Horvitz-Thompson estimator for the total of elements in the population.}
\item{VecStrataId.s}{vector of the strata identifiers; its length is equal to \eqn{n}, the sample size. Values in the argument \code{VecStrataId.s} must be numeric (integers) and must be sorted (increasingly). Note that all other arguments (variables, vectors, matrices) defined above that are of size \eqn{n} must follow the same order of \code{VecStrataId.s} correspondingly. Although this argument is not mandatory, if it is not provided the computations are made assuming that there is no stratification. There must not be any missing value.}
\item{VecStrataSizes.H}{vector of the strata population sizes; its length is equal to \eqn{H}, the number of strata. This vector contains the number of elements in the population for each stratum. The order of this vector must agree with the order of the argument \code{VecStrataId.s}. There must not be any missing value.}
\item{ShowStrata}{logical. If \code{TRUE} partial results from each stratum is displayed. This is an optional argument; default is \code{FALSE}.}
}
\details{
For the population total of the variable \eqn{y}:
\deqn{t = \sum_{k\in U} y_k}
the unbiased Narain (1951); Horvitz-Thompson (1952) estimator of \eqn{t} is given by:
\deqn{\hat{t}_{NHT} = \sum_{k\in s} w_k y_k}
where \eqn{w_k} denotes the sampling weight of the \eqn{k}-th element in the sample \eqn{s}, \eqn{w_k=1/\pi_k} with \eqn{\pi_k} denoting the inclusion probability of the \eqn{k}-th element in the sample. Let \eqn{\pi_{kl}} denotes the joint-inclusion probabilities of the \eqn{k}-th and \eqn{l}-th elements in the sample \eqn{s}. The variance of \eqn{\hat{t}_{NHT}} is given by:
\deqn{V(\hat{t}_{NHT}) = \sum_{k\in U}\sum_{l\in U} (\pi_{kl}-\pi_k\pi_l)w_k y_k w_l y_l}
which can therefore be estimated by the Horvitz-Thompson variance estimator (implemented by the current function if \code{VarEst="HT"}):
\deqn{\hat{V}(\hat{t}_{NHT}) = \sum_{k\in s}\sum_{l\in s} \frac{\pi_{kl}-\pi_k\pi_l}{\pi_{kl}}w_k y_k w_l y_l}

If the utilised sampling design is of fixed-size, the variance \eqn{V(\hat{t}_{NHT})} can be estimated by the Sen-Yates-Grundy variance estimator (implemented by the current function if \code{VarEst="SYG"}):
\deqn{\hat{V}(\hat{t}_{NHT}) = \frac{-1}{2}\sum_{k\in s}\sum_{l\in s} \frac{\pi_{kl}-\pi_k\pi_l}{\pi_{kl}}\left(w_k y_k - w_l y_l\right)^2}
For large-entropy sampling designs, the variance of \eqn{\hat{t}_{NHT}} is approximated by the Hajek (1964) variance:
\deqn{V(\hat{t}_{NHT}) \doteq \frac{N}{N-1}\left[\sum_{k\in U} w_k y_k^2\left(\frac{w_k-1}{w_k}\right)-dG^2\right]}
with \eqn{d=\sum_{k\in U}w_k^{-2}(w_k-1)} and \eqn{G=d^{-1}\sum_{k\in U}w_k^{-1}(w_k-1)y_k}.

This approximate variance can therefore be estimated by the variance estimator (implemented by the current function if \code{VarEst="Hajek"}):
\deqn{\hat{V}(\hat{t}_{NHT}) = \frac{n}{n-1}\left[\sum_{k\in s}w_k^2 y_k^2\left(\frac{w_k-1}{w_k}\right)-\hat{d}\hat{G}^2\right]}
where \eqn{\hat{d}=\sum_{k\in s}w_k^{-1}(w_k-1)} and \eqn{\hat{G}=\hat{d}^{-1}\sum_{k\in s}(w_k-1)y_k}.

The Hajek (1964) variance approximation is designed for large-entropy sampling designs and large populations, i.e. care should be taken with highly-stratified samples, e.g. Berger (2005).
  }
\value{
The function returns a dataframe with \eqn{Q} rows (the number of variables of interest) and some columns, whose number depends on the provided information and utilised expressions in computations. The results in the returned columns are:
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
\item{Nhat}{an estimate of the total of elements in population.}
\item{fhat}{an estimate fo the sampling fraction.}
\item{N}{the total of elements in population.}
\item{f}{the sampling fraction.}
If a stratified sampling design was specified and if \code{ShowStrata=TRUE} some further columns are displayed with partial results. Note that these per-stratum partial results are NOT returned by the function, they are only on-screen information.
\item{h}{stratum counter.}
\item{Stratum}{stratum identifier.}
\item{nh}{the sample size in the stratum \eqn{h}.}
\item{Nh}{the population total of elements in the stratum \eqn{h}.}
\item{fh}{the sampling fraction in the stratum \eqn{h}.}
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


##########################################
## Setting up sample data for the examples
##########################################

require(samplingVarEst)   #Loading the samplingVarEst package)
data(oaxaca)              #Loads the Oaxaca municipalities dataset (samplingVarEst)
MyPopulationData    <- oaxaca

#Data must be sorted according to stratification variable
MyPopulationData    <- MyPopulationData[order(MyPopulationData$IDREGION),]

#Defines the (already drawn) sample indicators to be used for the examples
SampleIndicators.U  <- MyPopulationData$sHOMES00        #Defines the sample indicators
N                   <- dim(MyPopulationData)[1]         #Defines the population size

#Reconstructs the population 1st order incl. probabilities (sample size was 373)
pik.U               <- Pk.PropNorm.U(373, MyPopulationData$HOMES00)
Weights.s           <- 1/pik.U[SampleIndicators.U==1]   #Creating the sample weights

#Approximates 2nd order incl. probs. using sample based quantitites
#Note that this approximation is only suitable for large-entropy sampling designs
pikl.s              <- Pkl.Hajek.s(pik.U[SampleIndicators.U==1])

#Defines the variables of interest y1 and y2 for those elements in the sample
y1.s                <- MyPopulationData$POP10[SampleIndicators.U==1]
y2.s                <- MyPopulationData$POPMAL10[SampleIndicators.U==1]

#Creates a matrix and a dataframe from the above observed variables of interest
MySampleMatrix      <- as.matrix(cbind(y1.s, y2.s), nrow= 373)
MySampleData        <- data.frame(cbind(y1.s, y2.s))

#Creates vectors of strata membership and strata sizes
StrataMembership    <- MyPopulationData$IDREGION[SampleIndicators.U==1]
StrataSizes         <- table(MyPopulationData$IDREGION)



#####################################################################
## Example 1: A variable of interest (vector), without stratification
#####################################################################

Estimate.Total.NHT(MatY.s  = y1.s     ,
                   VecWk.s = Weights.s)

Estimate.Total.NHT(MatY.s  = y1.s     ,
                   VecWk.s = Weights.s,
                   VarEst  = "HT"     )

Estimate.Total.NHT(MatY.s   = y1.s     ,
                   VecWk.s  = Weights.s,
                   VarEst   = "SYG"    ,
                   MatPkl.s = pikl.s   )

Estimate.Total.NHT(MatY.s   = y1.s     ,
                   VecWk.s  = Weights.s,
                   VarEst   = "SYG"    ,
                   MatPkl.s = pikl.s   ,
                   PopSize  = N        )

###################################################################################
## Example 2: A matrix/dataframe of 2 variables of interest, without stratification
###################################################################################

Estimate.Total.NHT(MatY.s   = MySampleMatrix,
                   VecWk.s  = Weights.s     ,
                   VarEst   = "SYG"         ,
                   MatPkl.s = pikl.s        ,
                   PopSize  = N             )

Estimate.Total.NHT(MatY.s   = MySampleData,
                   VecWk.s  = Weights.s   ,
                   VarEst   = "SYG"       ,
                   MatPkl.s = pikl.s      ,
                   PopSize  = N           )

##################################################################
## Example 3: A variable of interest (vector), with stratification
##################################################################

Estimate.Total.NHT(MatY.s           = y1.s            ,
                   VecWk.s          = Weights.s       ,
                   VecStrataId.s    = StrataMembership,
                   VecStrataSizes.H = StrataSizes     )

Estimate.Total.NHT(MatY.s           = y1.s            ,
                   VecWk.s          = Weights.s       ,
                   VecStrataId.s    = StrataMembership,
                   VecStrataSizes.H = StrataSizes     ,
                   ShowStrata       = TRUE            )

################################################################################
## Example 4: A matrix/dataframe of 2 variables of interest, with stratification
################################################################################

Estimate.Total.NHT(MatY.s           = MySampleData    ,
                   VecWk.s          = Weights.s       ,
                   VecStrataId.s    = StrataMembership,
                   VecStrataSizes.H = StrataSizes     ,
                   ShowStrata       = TRUE            )
}
\keyword{estimates}
\keyword{total}
