rm(list=ls())
options(digits=16)
# Example of NLS

# Data
#DATA=read.table("H12.csv",header=F,sep=',',col.names=c("t","F"))
t=c(1.1, 3.1, 5.1, 7.1, 12.1, 17.1, 22.1, 27.1, 32.1, 37.1, 42.1, 47.1, 52.1, 57.1, 62.1, 67.1, 72.1, 77.1, 82.1, 87.1, 92.1, 97.1, 102.1, 107.1, 115.1, 123.1, 131.1, 139.1, 147.1, 155.1, 163.1, 171.1, 179.1, 187.1, 195.1, 203.1, 211.1, 219.1, 227.1, 237.1, 247.1, 257.1, 267.1, 277.1, 287.1, 297.1, 307.1, 317.1, 327.1, 337.1, 347.1, 357.1, 367.1, 377.1, 387.1, 397.1, 407.1, 417.1, 427.1)
F=c(0.0977601998, 0.174783737, 0.2349572581, 0.2845963262, 0.3846295144, 0.4565562661, 0.5135314213, 0.5571324962, 0.5952466637, 0.6271657734, 0.6541027955, 0.6790957059, 0.7026772855, 0.7207867957, 0.7378093042, 0.7527797418, 0.7667195217, 0.7826460059, 0.7926127903, 0.8036890665, 0.8126971193, 0.8200648977, 0.8299734111, 0.8374219716, 0.8478726383, 0.8583356327, 0.8682122633, 0.8766409527, 0.8834150677, 0.8891601355, 0.8962832942, 0.9022124512, 0.9077612542, 0.914217411, 0.9207008732, 0.9224623019, 0.92808996, 0.9294081919, 0.9354150696, 0.9397991118, 0.9428710763, 0.9458128153, 0.9484635021, 0.9541724157, 0.9556748509, 0.9591829091, 0.9631354145, 0.9646526788, 0.9664148798, 0.9666498127, 0.9697618561, 0.9724616419, 0.973778977, 0.9750848225, 0.9792881657, 0.9793739011, 0.979256188, 0.979441914, 0.9797923127)
DATA=data.frame(cbind(t,F))

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

# How to eval data
RecoFunc2Pop( DATA[,1], c(0.071,0.00243,0.0184), c(0.719724956298514, 2.5, 100) )
RecoFunc3Pop( DATA[,1], c(25,0.0280,0.0072,0.3831,0.0541,0.0023,0.0009), c(0.719724956298514, 2.5, 100) )
# plot( DATA[,1], RecoFunc2Pop( DATA[,1], c(0.071,0.00243,0.0184), c(0.719724956298514, 2.5, 100) ) )
# plot( DATA[,1], RecoFunc3Pop( DATA[,1], c(25,0.0280,0.0072,0.3831,0.0541,0.0023,0.0009), c(0.719724956298514, 2.5, 100) ) )

# Fit model to data
fitM0D <- nls( F ~ RecoFunc2Pop(t, c(D, kb, ku), c(0.719724956298514, 2.5, 100) ), 
            data = DATA,
            start = c(D=0.071, kb = 0.0010819, ku = 0.0184), 
            trace = T,
            control = nls.control(warnOnly=TRUE,minFactor=0.0001, tol=1e-08),
            algorithm="port",
            lower = rep(0,2),
            upper = rep(1,2))
summary(fitM0D)

fitM1 <- nls( F ~ RecoFunc3Pop(t, c(25, kb, ku, gb, gu, eb, eu), c(0.719724956298514, 2.5, 100) ), 
            data = DATA,
            start = c(kb=0.0280, ku=0.0072, gb=0.3831, gu=0.0541, eb=0.0023,eu=0.0009), 
            trace = T,
            control = nls.control(warnOnly=TRUE,minFactor=0.000001, tol=1e-02),
            algorithm="port",
            lower = rep(0,6),
            upper = rep(1,6))
summary(fitM1)
print(coef(fitM1))

fitM9 <- nls( F ~ RecoFunc3Pop(t, c(25, 0, ku, gb, gu, eb, 0), c(0.719724956298514, 2.5, 100) ), 
              data = DATA,
              start = c(ku=0.00002, gb=0.00031, gu=0.0541, eb=0.0023), 
              trace = T,
              control = nls.control(warnOnly=TRUE,minFactor=0.000001, tol=1e-08),
              algorithm="port",
              lower = rep(0,6),
              upper = rep(1,6))
summary(fitM9)
print(coef(fitM9))
print(residuals(fitM9))
print(df.residual(fitM9))
print(deviance(fitM9))
print(logLik(fitM9))
