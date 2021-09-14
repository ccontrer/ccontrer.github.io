# Functions 
# All user-defined functions are here

from sage.rings.polynomial.real_roots import *
import numpy as np

def RecoFunc2Pop( t, par, theta ):

    r1 = theta[0]
    twoh = theta[1]
    N = theta[2]

    L = twoh / ( 1 - r1 ) 
    c = L / 2
    h = twoh / 2

    D = par[0]
    kappa_b = par[1]
    kappa_u = par[2]

    kappaub = kappa_b + kappa_u

    s = 0

    for n in range( 1, N ):
        DnpiL = D * ( n * float(pi) / L ) ** 2
        S = ( sin( n * float(pi) * ( c - h ) / L ) - sin( n * float( pi ) * ( c + h ) / L ) ) / ( n * float( pi ) )
        Q = sqrt( ( kappaub + DnpiL ) ^ 2 - 4 * kappa_u * DnpiL )
        r1 = - ( kappaub + DnpiL - Q ) / 2
        r2 = - ( kappaub + DnpiL + Q ) / 2
        F1 = ( r1 + kappaub ) / Q
        F2 = - ( r2 + kappaub ) / Q 
        B1 = kappaub + r1 + DnpiL
        B2 = kappaub + r2 + DnpiL
        
        s += ( B1 * F1 * exp( r1 * t ) + B2 * F2 * exp( r2 * t ) ) * S ** 2
    return 1 - ( L ^ 2 ) / ( h * ( L - 2 * h ) * ( kappaub ) ) * s  

####################################################################################################

def RecoFunc3Pop( t, par, theta ):
    
    r1 = theta[0]
    twoh = theta[1]
    N = theta[2]

    L = twoh / ( 1 - r1 ) 
    c = L / 2
    h = twoh / 2

    D = par[0]
    kappa_b = par[1]
    kappa_u = par[2]
    gamma_b = par[3]
    gamma_u = par[4]
    eta_b = par[5]
    eta_u = par[6]

    PI_u = kappa_u * eta_b + kappa_u * gamma_u + eta_u * gamma_u 
    PI_w = gamma_b * kappa_u + gamma_b * eta_u + kappa_b * eta_u 
    PI_v = eta_b * gamma_b + eta_b * kappa_b + gamma_u * kappa_b 
    PI_t = PI_u + PI_w + PI_v 

    Sigma_u = kappa_u + gamma_u + eta_b + eta_u 
    Sigma_t = kappa_b + kappa_u + gamma_b + gamma_u + eta_b + eta_u 
    
    rho = lambda ( s ): s ** 2 + Sigma_t * s + PI_t 

    F = 1

    for n in range( 1, N ):

         DnpiL = D * ( n * float( pi ) / L ) ** 2
         S = ( sin( n * float( pi ) * ( c - h ) / L ) - sin( n * float( pi ) * ( c + h ) / L ) ) / ( n * float( pi ) )
         z = Roots3deg( PI_u , PI_t , Sigma_u , Sigma_t , DnpiL ) #create this function
         r = map( lambda x: ( S * L ) ** 2 * rho( x ) / ( PI_t * ( L - 2 * h ) * h ), z )
         A = CoefLaplace( z , r )

         F -= ( ( PI_t + Sigma_u * DnpiL ) * LaplaceInv( t , z , A ) + ( Sigma_t + DnpiL ) * LaplaceInvD1( t , z , A ) + LaplaceInvD2( t , z , A) )
     
    return F

####################################################################################################

def LaplaceInv( t , z , A ):

    if not z[0].is_real():
        s = z[0].real
        w = z[0].imag
        a = abs( A[0] )
        c = arg( A[0] )
        L = 2 * a * exp( s * t ) * cos( w * t + c ) + A[2] * exp( z[2] * t )
    else:
        L = A[0] * exp( z[0] * t ) + A[1] * exp( z[1] * t ) + A[2] * exp( z[2] * t )
    return L

####################################################################################################

def LaplaceInvD1( t , z , A ):

    if not z[0].is_real():
        s = z[0].real
        w = z[0].imag
        a = abs( A[0] )
        c = arg( A[0] )
        L = 2 * a * exp( s * t ) * ( s * cos( w * t + c ) - w * sin( w * t + c ) ) + z[2] * A[2] * exp( z[2] * t )
    else:
        L = z[0] * A[0] * exp( z[0] * t ) + z[1] * A[1] * exp( z[1] * t ) + z[2] * A[2] * exp( z[2] * t )
    return L

####################################################################################################

def LaplaceInvD2( t , z , A ):

    if not z[0].is_real():
        s = z[0].real
        w = z[0].imag
        a = abs( A[0] )
        c = arg( A[0] )
        L = 2 * a * exp( s * t ) * ( ( s ** 2 - w ** 2 ) * cos( w * t + c ) - 2 * s * w * sin( w * t + c ) ) + z[2] ** 2 * A[2] * exp( z[2] * t )
    else:
        L = z[0] ** 2 * A[0] * exp( z[0] * t ) + z[1] ** 2 * A[1] * exp( z[1] * t ) + z[2] ** 2 * A[2] * exp( z[2] * t )
    return L

####################################################################################################

def CoefLaplace( z , r ):

    return [r[0] / ( ( z[0] - z[1] ) * ( z[0] - z[2] ) ),
            r[1] / ( ( z[1] - z[0] ) * ( z[1] - z[2] ) ), 
            r[2] / ( ( z[2] - z[0] ) * ( z[2] - z[1] ) )]
    

####################################################################################################

def Roots3deg( PI_u , PI_t , Sigma_u , Sigma_t , DnpiL ):

    #  Q = 1 / 9 * ( 3 * ( PI_t + Sigma_u * DnpiL ) - ( Sigma_t + DnpiL ) ^ 2 ) ;
    #  R = 1 / 54 * ( 9 * ( Sigma_t + DnpiL ) * ( PI_t + Sigma_u * DnpiL ) - 27 * PI_u * DnpiL - ( Sigma_t + DnpiL ) ^ 3 ) ;
    #  D = Q ^ 3 + R ^ 2 ;
    #  if D > 0
    #       S = sign( R + sqrt( D ) ) .* power( abs( R + sqrt( D ) ), 1 / 3 ) ;
    #       T = sign( R - sqrt( D ) ) .* power( abs( R - sqrt( D ) ), 1 / 3 ) ;
    #  elseif D <= 0
    #       rho = sqrt( - ( Q ^ 3 ) ) ;
    #       phi = acos( R / rho );
    #       S = sign( rho ) .* power( abs( rho ), 1 / 3 ) * cos( phi / 3 ) ;
    #       T = sign( rho ) .* power( abs( rho ), 1 / 3 ) * cos( - phi / 3 ) ; ;
    #  end
    #  z( 1 ) = - 0.5 * (S + T) - 1 / 3 * (Sigma_t + DnpiL) + 0.5 * i * sqrt( 3 ) * ( S - T ) ;
    #  z( 2 ) = - 0.5 * (S + T) - 1 / 3 * (Sigma_t + DnpiL) - 0.5 * i * sqrt( 3 ) * ( S - T ) ;
    #  z( 3 ) = S + T - 1 / 3 * ( Sigma_t + DnpiL ) ;

    #t = polygen( RR )
    #z = ( t ** 3 + ( Sigma_t + DnpiL ) * t ** 2 + ( PI_t + Sigma_u * DnpiL ) * t + ( PI_u * DnpiL ) ).complex_roots()
    z = list( np.roots( [ 1, Sigma_t + DnpiL, PI_t + Sigma_u * DnpiL, PI_u * DnpiL ] ) )
    z[0] = RR( z[0] )
    z[1] = RR( z[1] )
    z[2] = RR( z[2] )
    if z[2].imag > 1e-14:
         temp = z[0]
         z[0] = z[2]
         z[2] = temp
    return z

