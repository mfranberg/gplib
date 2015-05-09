\name{predict.gp}
\alias{predict.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.gp(object, newdata, level = 0.95, interval = c("none", "prediction"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
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
function (object, newdata, level = 0.95, interval = c("none", 
    "prediction")) 
{
    if (missing(newdata)) {
        Xstar = model.matrix(object)
    }
    else {
        tt = terms(object$model)
        tr = delete.response(tt)
        trd = drop.terms(tr)
        mf = model.frame(trd, newdata)
        Xstar = model.matrix(trd, mf)
    }
    X = model.matrix(object)
    Kstar = matrix(0, nrow(Xstar), nrow(X))
    for (i in 1:nrow(Xstar)) {
        for (j in 1:nrow(X)) {
            Kstar[i, j] = object$kernel$kernelfunc(Xstar[i, ], 
                X[j, ], object$kernel$params)
        }
    }
    fbar = Kstar \%*\% object$alpha
    V = rep(0, nrow(Xstar))
    for (i in 1:nrow(Xstar)) {
        v = solve(t(object$L), Kstar[i, ])
        kstar = object$kernel$kernelfunc(Xstar[i, ], Xstar[i, 
            ], object$kernel$params)
        V[i] = kstar - t(v) \%*\% v
    }
    interval = match.arg(interval)
    if (interval == "none") {
        return(fbar)
    }
    else if (interval == "prediction") {
        lower = qnorm(1 - level, fbar, V)
        upper = qnorm(level, fbar, V)
        ret = cbind(fbar, lower, upper)
        colnames(ret) = c("fit", "lwr", "upr")
        return(ret)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
