# Function to evaluate (model) usage recoFunM0( dataTime, parametersToFit, otherParameters )
RecoFunc2Pop <- function( t, par, theta )
{  
    r1 = theta[1]
    twoh = theta[2]
    N = theta[3]
    
    L = twoh / ( 1 - r1 ) 
    c = L / 2
    h = twoh / 2
    
    D = par[1]
    kappa_b = par[2]
    kappa_u = par[3]
    
    kappaub = kappa_b + kappa_u
    
    s = 0
    
    for( n in 1:N ){
      DnpiL = D * ( n * pi / L ) ** 2
      S = ( sin( n * pi * ( c - h ) / L ) - sin( n * pi * ( c + h ) / L ) ) / ( n * pi )
      Q = sqrt( ( kappaub + DnpiL ) ^ 2 - 4 * kappa_u * DnpiL )
      r1 = - ( kappaub + DnpiL - Q ) / 2
      r2 = - ( kappaub + DnpiL + Q ) / 2
      F1 = ( r1 + kappaub ) / Q
      F2 = - ( r2 + kappaub ) / Q 
      B1 = kappaub + r1 + DnpiL
      B2 = kappaub + r2 + DnpiL
      s = s + ( B1 * F1 * exp( r1 * t ) + B2 * F2 * exp( r2 * t ) ) * S ^ 2
    }
    
    return(1 - ( L ^ 2 ) / ( h * ( L - 2 * h ) * ( kappaub ) ) * s)
}

RecoFunc3Pop <- function( t , par , theta )
{
    r1 = theta[1]
    twoh = theta[2]
    N = theta[3]
    
    L <- twoh / ( 1 - r1 ) 
    c <- L / 2 
    h <- twoh / 2 
    
    D <- par[ 1 ] 
    kappa_b <- par[ 2 ] 
    kappa_u <- par[ 3 ] 
    gamma_b <- par[ 4 ] 
    gamma_u <- par[ 5 ] 
    eta_b <- par[ 6 ] 
    eta_u <- par[ 7 ]  
        
    PI_u <- kappa_u * eta_b + kappa_u * gamma_u + eta_u * gamma_u 
    PI_w <- gamma_b * kappa_u + gamma_b * eta_u + kappa_b * eta_u 
    PI_v <- eta_b * gamma_b + eta_b * kappa_b + gamma_u * kappa_b 
    PI_t <- PI_u + PI_w + PI_v 
    Sigma_u <- kappa_u + gamma_u + eta_b + eta_u 
    Sigma_t <- kappa_b + kappa_u + gamma_b + gamma_u + eta_b + eta_u 

    rho <- function( s ) s ^ 2 + Sigma_t * s + PI_t
  
    F <- 1 
    for ( n in 1 : N )
    {
        DnpiL <- D * ( n * pi / L ) ^ 2 
        S <- ( sin( n * pi * ( c - h ) / L ) - sin( n * pi * ( c + h ) / L ) ) / ( n * pi ) 
        z <- Roots3deg( PI_u , PI_t , Sigma_u , Sigma_t , DnpiL )  
        r <- ( S * L ) ^ 2 * rho( z ) / ( PI_t * ( L - 2 * h) * h ) 
        A <- CoefLaplace( z , r ) 
        F <- F - ( ( PI_t + Sigma_u * DnpiL ) * LaplaceInv( t , z , A ) + 
                       ( Sigma_t + DnpiL ) * LaplaceInvD1( t , z , A ) + LaplaceInvD2( t , z , A) ) 
    }
    return( F )
}

LaplaceInv <- function( t , z , A )
{
    if (! is.real( z[ 1 ] ) )
    {
        s = Re( z[ 1 ] ) 
        w = Im( z[ 1 ] ) 
        a = abs( A[ 1 ] ) 
        c = Arg( A[ 1 ] ) 
        L = 2 * a * exp( s * t ) * cos( w * t + c ) + A[ 3 ] * exp( z[ 3 ] * t ) 
    }
    else
    {
        L = A[ 1 ] * exp( z[ 1 ] * t ) + A[ 2 ] * exp( z[ 2 ] * t ) + A[ 3 ] * exp( z[ 3 ] * t ) 
    }
    return( L )
}

LaplaceInvD1 <- function( t , z , A )
{
    if (! is.real( z[ 1 ] ) )
    {
        s = Re( z[ 1 ] ) 
        w = Im( z[ 1 ] ) 
        a = abs( A[ 1 ] ) 
        c = Arg( A[ 1 ] ) 
        L = 2 * a * exp( s * t ) * ( s * cos( w * t + c ) - w * sin( w * t + c ) ) + z( 3 ) * A( 3 ) * exp( z( 3 ) * t ) 
    }
    else
    {
        L <- z[ 1 ] * A[ 1 ] * exp( z[ 1 ] * t ) + z[ 2 ] * A[ 2 ] * exp( z[ 2 ] * t ) + z[ 3 ] * A[ 3 ] * exp( z[ 3 ] * t )
    }
    return( L )
}

LaplaceInvD2 <- function( t , z , A )
{
    if (! is.real( z[ 1 ] ) )
    {
        s = Re( z[ 1 ] ) 
        w = Im( z[ 1 ] ) 
        a = abs( A[ 1 ] ) 
        c = Arg( A[ 1 ] ) 
        L = 2 * a * exp( s * t ) * ( ( s ^ 2 - w ^ 2 ) * cos( w * t + c ) - 2 * s * w * sin( w * t + c ) ) +
            z( 3 ) ^ 2 * A( 3 ) * exp( z( 3 ) * t ) 
    }
    else
    {
        L <- z[ 1 ] ^ 2 * A[ 1 ] * exp( z[ 1 ] * t ) + z[ 2 ] ^ 2 * A[ 2 ] * exp( z[ 2 ] * t ) +
            z[ 3 ] ^ 2 * A[ 3 ] * exp( z[ 3 ] * t )
    }
    return( L )
}

CoefLaplace <- function( z , r )
{
    A1 <- r[ 1 ] / ( ( z[ 1 ] - z[ 2 ] ) * ( z[ 1 ] - z[ 3 ] ) )   
    A2 <- r[ 2 ] / ( ( z[ 2 ] - z[ 1 ] ) * ( z[ 2 ] - z[ 3 ] ) )   
    A3 <- r[ 3 ] / ( ( z[ 3 ] - z[ 1 ] ) * ( z[ 3 ] - z[ 2 ] ) ) 
    return( c( A1 , A2 , A3 ) )
}

Roots3deg <- function( PI_u , PI_t , Sigma_u , Sigma_t , DnpiL )
{
    z <- polyroot( c( PI_u * DnpiL ,  PI_t + Sigma_u * DnpiL , Sigma_t + DnpiL , 1 ) )
    if ( Im( z[ 3 ] ) > 1e-14 )
    {
        t = z[ 1 ]
        z[ 1 ] = z[ 3 ]
        z[ 3 ] = t
    }
    else{
        z = Re(z)
    }
    return( z )
}