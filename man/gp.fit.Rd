\name{gp.fit}
\alias{gp.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gp.fit(X, y, kernel = kernel, sigmastart = 1, optimize = True)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
}
  \item{sigmastart}{
%%     ~~Describe \code{sigmastart} here~~
}
  \item{optimize}{
%%     ~~Describe \code{optimize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, kernel = kernel, sigmastart = 1, optimize = True) 
{
    optimize_gp = function(X, y, kernel, kernel.params, sigma) {
        n = nrow(X)
        K = matrix(0, ncol = n, nrow = n)
        for (i in 1:n) {
            for (j in 1:n) {
                K[i, j] = kernel$kernelfunc(X[i, ], X[j, ], kernel.params)
            }
        }
        L = chol(K + diag(sigma^2, n))
        alpha = solve(L, solve(t(L), y))
        loglik = -0.5 * t(y) \%*\% alpha - sum(diag(L)) - 0.5 * 
            n * log(2 * pi)
        list(L = L, alpha = alpha, loglik = loglik, sigma = sigma)
    }
    gradient_gp = function(p, Kinv, alpha, kernel) {
        0.5 * tr(alpha \%*\% t(alpha) - Kinv) * kernel$gradient(p)
    }
    nlm_gp = function(p, X, y, kernel) {
        -optimize_gp(X, y, kernel, p[2:(length(kernel$params) + 
            1)], p[1])$loglik
    }
    params = nlm(nlm_gp, c(1, kernel$params), X = X, y = y, kernel = kernel)
    kernel$params = params[2:length(params)]
    optimize_gp(X, y, kernel, kernel$params, params[1])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
