make.kernel = function(kernel)
{
    switch(kernel, squared_exponential = {
        name = kernel
        params = c( 1.0, 1.0 )
        params.names = c( "sigmak", "lengthscale" )
        kernelfunc = function(x1, x2, kernel.params) kernel.params[ 1 ] * exp( -0.5*sum( (x1-x2)^2 ) / kernel.params[ 2 ] )
    }, stop( gettextf( "%s kernel not recognised", sQuote( kernel ) ), domain = NA ) )
    
    structure( list( name = name,
                     params = params,
                     params.names = params.names,
                     kernelfunc = kernelfunc ), class = "kernel-gp" )
}

gp = function(formula, data, kernel = make.kernel( "squared_exponential" ) )
{
    if( missing( data ) )
    {
        data = environment( formula )
    }
    
    m = model.frame( formula, data )
    X = model.matrix( formula, m )
    y = model.response( m )
    
    fit = gp.fit( X, y, kernel )
    
    z = c( fit, list( model = m ) )
    class( z ) = "gp"
    
    return( z )
}

gp.fit = function(X, y, kernel = kernel, sigmastart = 1.0, optimize = True)
{   
    optimize_gp = function(X, y, kernel, kernel.params, sigma)
    {
        n = nrow( X )
        K = matrix( 0, ncol = n, nrow = n )
        for(i in 1:n)
        {
            for(j in 1:n)
            {
                K[ i, j ] = kernel$kernelfunc( X[ i, ], X[ j, ], kernel.params )
            }
        }
        
        L = chol( K + diag( sigma^2, n ) )
        alpha = solve( L, solve( t( L ), y ) )
        loglik = -0.5 * t( y ) %*% alpha - sum( diag( L ) ) - 0.5 * n * log( 2 * pi )
        
        list( L = L, alpha = alpha, loglik = loglik, sigma = sigma )
    }
    
    gradient_gp = function(p, Kinv, alpha, kernel)
    {
        0.5 * tr( alpha %*% t( alpha ) - Kinv ) * kernel$gradient( p )
    }
    
    nlm_gp = function( p, X, y, kernel )
    {
        -optimize_gp( X, y, kernel, p[ 2:(length( kernel$params )+1) ], p[ 1 ] )$loglik
    }
    
    params = nlm( nlm_gp, c( 1.0, kernel$params ), X = X, y = y, kernel = kernel )
    kernel$params = params[ 2:length( params ) ]
    
    optimize_gp( X, y, kernel, kernel$params, params[ 1 ] )
}

model.matrix.gp = function(object)
{
    model.matrix( object$model )
}

print.gp = function(x)
{
    cat( "Call:\n" )
    print( )
    cat( "Marginal likelihood:", x$loglik, "\n" )
    cat( "sigma:", x$sigma, "\n" )
    kernel.params = x$kernel$params
    names( kernel.params ) = x$kernel$params.names
    
    cat( "kernel:", x$kernel$name , "\n" )
    cat( "kernel params:\n" )
    print( kernel.params )
    invisible( x )
}

predict.gp = function(object, newdata, level = 0.95, interval = c( "none", "prediction" ) )
{
    if( missing( newdata ) )
    {
        Xstar = model.matrix( object )
    } else
    {
        tt = terms( object$model )
        tr = delete.response( tt )
        trd = drop.terms( tr )
        mf = model.frame( trd, newdata )
        Xstar = model.matrix( trd, mf )
    }
    X = model.matrix( object )
    Kstar = matrix( 0, nrow( Xstar ), nrow( X ) )
    for(i in 1:nrow( Xstar ) )
    {
        for(j in 1:nrow( X ) )
        {
            Kstar[ i, j ] = object$kernel$kernelfunc( Xstar[ i, ], X[ j, ], object$kernel$params )
        }
    }
    
    fbar = Kstar %*% object$alpha
    V = rep( 0, nrow( Xstar ) )
    for(i in 1:nrow( Xstar ) )
    {
        v = solve( t( object$L ), Kstar[ i, ] )
        kstar = object$kernel$kernelfunc( Xstar[ i, ], Xstar[ i, ], object$kernel$params )
        V[ i ] = kstar - t( v ) %*% v
    }
    
    interval = match.arg( interval )
    if( interval == "none" )
    {
        return( fbar )
    }
    else if( interval == "prediction" )
    {
        lower = qnorm( 1 - level, fbar, V )
        upper = qnorm( level, fbar, V )
        
        ret = cbind( fbar, lower, upper )
        colnames( ret ) = c( "fit", "lwr", "upr" )
        
        return( ret )
    }
}

#x = rnorm( 100, 0, 5 )
#y = sin( 0.5*x )
#ynoise = sin( 0.5*x ) + rnorm( 100, 0, 0.1)
#original = data.frame( x = x, y = y, ynoise = ynoise )

#stats = gp( ynoise ~ x, data = original )

#xnew = seq( min( x ), max( x ), 0.5 )
#ynew = predict( stats, data.frame( x = xnew ) )

#prediction = as.data.frame( predict( stats, data.frame( x = xnew ), interval = "prediction" ) )
#prediction$x = xnew

#library( ggplot2 )

#ggplot( original, aes( x = x, y = y ) ) + geom_line( ) +
#    geom_point( aes( x = x, y = ynoise ) ) +
#    geom_line( data = prediction, aes( x = x, y = fit ), color = "blue" ) +
#    geom_ribbon( data = prediction, aes( ymin = lwr, ymax = upr, x = x, y = fit ), alpha = 0.2 )

