\name{make.kernel}
\alias{make.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.kernel(kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
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
function (kernel) 
{
    switch(kernel, squared_exponential = {
        name = kernel
        params = c(1, 1)
        params.names = c("sigmak", "lengthscale")
        kernelfunc = function(x1, x2, kernel.params) kernel.params[1] * 
            exp(-0.5 * sum((x1 - x2)^2)/kernel.params[2])
    }, stop(gettextf("\%s kernel not recognised", sQuote(kernel)), 
        domain = NA))
    structure(list(name = name, params = params, params.names = params.names, 
        kernelfunc = kernelfunc), class = "kernel-gp")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
