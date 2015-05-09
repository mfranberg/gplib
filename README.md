# gplib

The gplib library is a simplistic implementation of Guassian Process regression in R. It
can be used in a similar way to other models in R, for example

    library( gplib )
    
    gp( y ~ x1 + x2 + x3, data = data, kernel = make.kernel( "squared_exponential" ) )

And the return value supports some of the normal accessor functions like residuals, predict
and print.
