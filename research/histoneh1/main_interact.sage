load_attach_path('/home/shigoto/Documents/Studies/PhD/Courses/MATH 600 Sage course (computing in mathematics)/Test/')
load "data.sage"
load "functions.sage"

# GUI
import StringIO
import csv
@interact
def recovery(
             data=input_box( default=Data, width=20, height=10, type=str ),
             kappa_b = slider( srange(0, 0.2, 0.00001), default = 0.00243, label='$\kappa_b$' ),
             kappa_u = slider( srange(0, 0.2, 0.00001), default = 0.0184, label='$\kappa_u$' ),
             D = slider( srange(0, 0.1, 0.0001), default = 0.071, label='$D$' ),
             fit_bool = checkbox( False, label='Fit model to data with initial parameters', )
             ):
             
    input_file = StringIO.StringIO(data.strip())    
    data = csv.reader(input_file)
    x = [], y = []
    for row in data:
        try:
           temp=map(RDF, row)
           x.append(temp[0])
           y.append(temp[1])
        except ValueError:
            continue
    
    gr = point( zip(x, y), ymin=0, ymax=1, color='green', legend_label='Data' )
    t = var('t')
    if fit_bool == False:
        gr += plot( recoFunM0( t, [kappa_b, kappa_u, D], [0.719724956298514, 2.5, 100]), (t, 0, max(x)), legend_label='Fit (initial parameters)' ) 
    if fit_bool == True:
        var( 't, kappa_b_hat, kappa_u_hat, D_hat' )
        p_hat = find_fit(zip(x, y), 
                         recoFunM0(t, [kappa_b_hat, kappa_u_hat, D_hat], [0.719724956298514, 2.5, 100]), 
                         initial_guess = [kappa_b, kappa_u, D], 
                         parameters = [kappa_b_hat, kappa_u_hat, D_hat], 
                         variables = [t], 
                         solution_dict = True)
        gr += plot( recoFunM0(t, [p_hat[kappa_b_hat], p_hat[kappa_u_hat], p_hat[D_hat]], [0.719724956298514, 2.5, 100]), (t, 0, max(x)), legend_label='Fit' )  
         
    gr.axes_labels( ['$t$ time (sec.)', '$F(t)$ Fluorescence recovery'] )
    gr.set_legend_options(loc=4, borderaxespad=2 )
    gr.show()