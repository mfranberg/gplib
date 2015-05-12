library( ggplot2 )

make.kernel = function(kernel, kernel.params)
{
    switch(kernel, squared_exponential = {
        name = kernel
        params = c( 0, 0, 0 )
        if( !missing( kernel.params ) )
        {
            params = log( kernel.params )
        }
        
        params.names = c( "sigma", "sigmak", "lengthscale" )
        params.format = function( p ) c( sqrt( exp( p[ 1 ] ) ), exp( p[ 2:3 ] ) )
        kernelfunc = function(x1, x2, kernel.params, same = FALSE) exp( kernel.params[ 2 ] ) * exp( -0.5*sum( (x1-x2)^2 ) / exp( kernel.params[ 3 ] ) ) + same * exp( kernel.params[ 1 ] )
        kernelgrad = function(x1, x2, kernel.params, same = FALSE) c( same * exp( kernel.params[ 1 ] ),
                                                                      2*exp( kernel.params[ 2 ] ) * exp( -0.5*sum( (x1-x2)^2 ) / exp( kernel.params[ 3 ] )^2 ),
                                                                      exp( kernel.params[ 2 ] ) * exp( -0.5*sum( (x1-x2)^2 ) / exp( kernel.params[ 3 ] )^2 ) * -0.5*sum( (x1-x2)^2 ) / exp( kernel.params[ 3 ] )^2 )
    }, stop( gettextf( "%s kernel not recognised", sQuote( kernel ) ), domain = NA ) )
    
    structure( list( name = name,
                     params = params,
                     params.names = params.names,
                     params.format = params.format,
                     kernelfunc = kernelfunc,
                     kernelgrad = kernelgrad ), class = "kernel-gp" )
}

gp = function(formula, data, kernel = make.kernel( "squared_exponential" ), optimize = TRUE, startp = kernel$params )
{
    if( missing( data ) )
    {
        data = environment( formula )
    }
    
    mf = match.call( expand.dots = FALSE )
    m = match( c( "formula", "data" ), names( mf ), 0L )
    mf = mf[ c( 1L, m ) ]
    mf[[ 1L ]] = quote( stats::model.frame )
    mf = eval( mf, parent.frame( ) )
    
    X = model.matrix( formula, mf )
    y = model.response( mf )
    
    fit = gp.fit( X, y, kernel, optimize = optimize, startp = startp )
    
    z = c( fit, list( model = mf ) )
    class( z ) = "gp"
    
    return( z )
}

gp.compute_K = function(X, kernel, kernel.params)
{
    n = nrow( X )
    K = matrix( 0, ncol = n, nrow = n )
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            K[ i, j ] = kernel$kernelfunc( X[ i, ], X[ j, ], kernel.params, same = i == j )
        }
    }
    return( K )
}

gp.compute_Kgrad = function(X, kernel, kernel.params)
{
    n = nrow( X )
    Klist = list( )
    for(i in 1:length(kernel.params))
    {
        Klist[[ i ]] = matrix( 0, ncol = n, nrow = n )
    }
    
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            kgrad = -kernel$kernelgrad( X[ i, ], X[ j, ], kernel.params, i == j )
            for(k in 1:length(kernel.params))
            {
                Klist[[ k ]][ i, j ] = kgrad[ k ]
            }
        }
    }
    
    return( Klist )
}

gp.fit = function(X, y, kernel = kernel, optimize = TRUE, startp = kernel$params )
{   
    optimize_gp = function(X, y, kernel, kernel.params)
    {
        n = nrow( X )
        K = gp.compute_K( X, kernel, kernel.params )
        
        L = chol( K )
        alpha = solve( L, solve( t( L ), y ) )
        loglik = -0.5 * t( y ) %*% alpha - sum( diag( L ) ) - 0.5 * n * log( 2 * pi )
        
        list( L = L, alpha = alpha, loglik = loglik, kernel = kernel )
    }
    
    gradient_gp = function(p, X, y, kernel)
    {
        grad = rep( 0, length( p ) )
        
        K = gp.compute_K( X, kernel, p )
        Kinv = solve( K )
        Kgrad = gp.compute_Kgrad( X, kernel, p )
        
        alpha = Kinv %*% y
        for(i in 1:length(p))
        {
            grad[ i ] = 0.5 * sum( diag( ( alpha %*% t( alpha ) - Kinv ) %*% Kgrad[[ i ]] ) )
        }
        
        return( grad )
    }
    
    f_gp = function(p, X, y, kernel)
    {
        -optimize_gp( X, y, kernel, p )$loglik
    }
    
    if( optimize )
    {
        params = optim( startp, f_gp, gradient_gp, method = "CG", X = X, y = y, kernel = kernel, control = list( trace = 1, maxit = 5, reltol = 1e-1 ) )
        kernel$params = params$par
    }
    
    optimize_gp( X, y, kernel, kernel$params )
}

model.matrix.gp = function(object)
{
    model.matrix( object$model )
}

print.gp = function(x)
{
    cat( "Marginal likelihood:", x$loglik, "\n" )
    kernel.params = x$kernel$params
    names( kernel.params ) = x$kernel$params.names
    
    cat( "kernel:", x$kernel$name , "\n" )
    cat( "kernel params:\n" )
    print( x$kernel$params.format( kernel.params ) )
    invisible( x )
}

logLik.gp = function(object)
{
    object$loglik
}

predict.gp = function(object, newdata, X, level = 0.95, interval = c( "none", "prediction" ) )
{
    if( !missing( newdata ) )
    {
        tt = terms( object$model )
        tr = delete.response( tt )
        trd = drop.terms( tr )
        mf = model.frame( trd, newdata )
        Xstar = model.matrix( trd, mf )
    } 
    else if( !missing( X ) )
    {
        Xstar = X
    }
    else
    {
        Xstar = model.matrix( object )
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
        lower = qnorm( 1 - level, fbar, sqrt( V ) )
        upper = qnorm( level, fbar, sqrt( V ) )
        
        ret = cbind( fbar, lower, upper )
        colnames( ret ) = c( "fit", "lwr", "upr" )
        
        return( ret )
    }
}

plot.gp = function(object, true_y = NULL)
{
    tt = terms( object$model )
    X = model.matrix( tt )
    yold = model.response( tt )
    if( ncol( X ) == 1 || ( ncol( X ) == 2 && all( X[ ,1 ] == 1 ) ) )
    {
        xname = colnames( X )[ 1 ]
        xold = X[, 1 ]
        if( ncol( X ) == 2 )
        {
            xname = colnames( X )[ 2 ]
            xold = X[ ,2 ]
        }
        
        xnew = seq( min( xold ), max( xold ), (max( xold ) - min( xold ))/1000 )
        new_data = data.frame( x = xnew )
        colnames( new_data ) = xname
        
        old_data = data.frame( x = xold, y = yold )
        
        prediction = as.data.frame( predict( object, new_data, interval = "prediction" ) )
        prediction$x = xnew

        p = ggplot( old_data, aes( x = x, y = y ) ) +
            geom_point( ) +
            geom_line( data = prediction, aes( x = x, y = fit ), color = "blue" ) +
            geom_ribbon( data = prediction, aes( ymin = lwr, ymax = upr, x = x, y = fit ), alpha = 0.2 ) +
            scale_x_continuous( xname ) + 
            scale_y_continuous( colnames( yold ) )
        if( !missing( true_y ) )
        {
            prediction$ytrue = true_y( xnew )
            p = p + geom_line( data = prediction, aes( x = x, y = ytrue ) ) 
        }
        p
    }
    else
    {
        prediction = as.data.frame( predict( object, interval = "prediction" ) )
        prediction$y = yold
        prediction$residual = prediction$fit - yold
        ggplot( prediction, aes( x = fit, y = residual ) ) +
            geom_point(  ) + geom_abline( intercept = 0, slope = 0 )
    }
}

test = function()
{
    x = -10 + 20 * runif( 100 )
    y = sin( 0.5*x ) + rnorm( 100, 0, 0.5 )
    stats = gp( y ~ x )
    plot( stats, true_y = function(x) sin(0.5*x) )
}