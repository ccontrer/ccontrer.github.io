load_attach_path('/home/shigoto/Documents/Studies/PhD/Courses/MATH 600 Sage course (computing in mathematics)/Test/')
load "data.sage"
load "functions.sage"

kappa_b = 0.00243
kappa_u = 0.0184
D = 0.071
var( 't, kappa_b_hat, kappa_u_hat, D_hat' )
p_hat = find_fit(zip(x, y), 
                 recoFunM0(t, [kappa_b_hat, kappa_u_hat, D_hat], [0.719724956298514, 2.5, 100]), 
                 initial_guess = [kappa_b, kappa_u, D], 
                 parameters = [kappa_b_hat, kappa_u_hat, D_hat], 
                 variables = [t], 
                 solution_dict = True)
print "\n\n"
print p_hat
