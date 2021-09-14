load_attach_path('/home/shigoto/Documents/Studies/PhD/Courses/MATH 600 Sage course (computing in mathematics)/Test/')
%runfile "data.sage"
%attach "functions.sage"



RecoFunc2Pop(x[0], [0.071, 0.00243, 0.0184], [0.719724956298514, 2.5, 100])
RecoFunc3Pop(x[0], [0.071, 0.00243, 0.0184, 0.00243, 0.00346, 0.003463, 0.003526], [0.719724956298514, 2.5, 100])


#kappa_b = 0.00243
#kappa_u = 0.0184
#D = 0.071
#var( 't, kappa_b_hat, kappa_u_hat, D_hat' )
#p_hat = find_fit(zip(x, y), 
                 #recoFunM0(t, [kappa_b_hat, kappa_u_hat, D_hat], [0.719724956298514, 2.5, 100]), 
                 #initial_guess = [kappa_b, kappa_u, D], 
                 #parameters = [kappa_b_hat, kappa_u_hat, D_hat], 
                 #variables = [t], 
                 #solution_dict = True)
#print "\n\n"
#print p_hat
#print "\n"

#y_t = []
#for i in range(len(x)):
    #y_t.append( recoFunM0(x[i], [0.00243, 0.0184, 0.071], [0.719724956298514, 2.5, 100]) )

#from scipy.optimize import minimize, rosen, rosen_der

#def mintofit( fun, x, y, p, theta):
    #res = 0
    #for i in range(len(x)):
        #res += ( y[i] - fun(x[i], p, theta) ) ** 2
    #return res

#print mintofit( recoFunM0 , x , y, [0.00243, 0.0184, 0.071], [0.719724956298514, 2.5, 100] )

#attemp to use minimize as a fitting method
# the function mintofit computes the residial for |y-f(x;par)| so its minimizable
#bnds = ((0, 1), (0, 1))
#fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
#fun1 = lambda x,a: (a[0]*x[0] - 1)**2 + (a[1]*x[1] - 2.5)**2
#fun2 = lambda x: fun1(x,[3,2])
#res = minimize(lambda x: fun1(x,[3,2]), (2, 0), method='SLSQP', bounds=bnds)
#print res
#res = minimize( lambda p: mintofit( recoFunM0, x, y, p, [0.719724956298514, 2.5, 100] ), 
               #[0.00243, 0.0184, 0.071], 
               #method = 'TNC', 
               #bounds = ((0, 1), (0, 1), (0, 1)),
               #options={'maxiter': 100} )
#print res