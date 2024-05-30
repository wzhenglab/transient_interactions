import numpy as np                                                              
import os,sys                                                                   
import scipy.optimize as op                                                     
import matplotlib.pyplot as plt

                                                                                
"""                                                                             
fit data with adjusted polymer model and return ft nu                               

make sure directory contains polymer.py and a file of 
    labelling distances per frame

usage
    python3 fit_pr.py dye.npy
"""                                                                             
    

# polymer module contains all equations for the adjusted polymer model
import polymer as pol                                                           
def fit_func(params,prs,rc):                                                    
    """
    adjusted p(r) distribution has four free parameters
        ft: fraction of transient interactions
        nu: polymer scaling exponent
        xc: position of short range peak
        sigma: width of short range peak
    loss fuction to be minimized by fit_pr function
    """
    w = params[0]                                                               
    nu = params[1]                                                              
    xc = params[2]                                                              
    sigma = params[3]                                                           
    b = (0.55*99**nu)**2                                                        
    # calculate chi squared function between adjusted model and real p(r)
    guess = pol.pr_saw_gau(rc,w,b,nu,xc,sigma)                                  
    chi=np.sum((prs-guess)**2)                                                  
    return chi                                                                  
 

def fit_pr(prs,rc):
    """
    arguments
        p(r) to be fitted
        bins of p(r) distribution
    outputs 
        array of fitted parameters from scipy.optimize.minimize
    """
    adds = (prs,rc)     
    # define intial guesses
    # assume 0 transient interactions, excluded volume scaling, and positions/width
    #   guess come from behavior of model polymer electrostatic interactions
    x0=[0,0.553,1,0.6]
    # SLSQP seems to be most reliable fitting method, but feel free to optimize
    minum = op.minimize(fit_func,x0=x0,args=adds,#bounds=bounds,                 
            options={'disp':True},method='SLSQP')                               
    return minum.x                                                              


if __name__ == "__main__": 
    # bins used for making p(r) distribution
    # distance units are in nm
    # width of p(r) bins
    dx=0.1
    # 20 nm is generally large enough to accomodate labelling positions <200 residues apart
    x=np.arange(0,20,dx)
    xc=(x[1:]+x[:-1])/2
    # load the array of labelling ditsances per frame
    dye_file=sys.argv[1]
    dye=np.load(dye_file)
    # calculate p(r) distribution from labelling distances
    pr,q=np.histogram(dye,density=True,bins=x)
    #rescale p(r) so that sum=1
    pr=pr*dx
    fit_params=fit_pr(pr,xc)

var=['f_t','nu','xc','sigma']
for i,v in enumerate(var):
    print(v,'=',str(fit_params[i]))

