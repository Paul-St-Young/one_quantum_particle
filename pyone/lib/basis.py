import numpy as np
from scipy.special import factorial

def gto_norm(l,m,n,alpha):
    return (2*alpha/np.pi)**(3./4) * (\
(8.*alpha)**(l+m+n) *factorial(l) *factorial(m) *factorial(n)/\
(factorial(2*l) *factorial(2*m) *factorial(2*n))\
)**(1./2) 
# end def gto_norm

def gto_dict(x,y,z,params):
    # params should be a dictionary of parameters
    l = params["l"]
    m = params["m"]
    n = params["n"]
    alpha = params["alpha"]
    xo    = params["xo"]
    yo    = params["yo"]
    zo    = params["zo"]
    
    norm = gto_norm(l,m,n,alpha)
    r2 = (x-xo)*(x-xo) + (y-yo)*(y-yo) + (z-zo)*(z-zo)
    return norm *(x-xo)**l *(y-yo)**m *(z-zo)**n *np.exp(-alpha*r2)
# end def gto_dict

def et20_basis():
    """ Modified ET-555 basis set. Removed the most diffuse functions with alpha=2.091532.
    Only include basis functions that are symmetric around the z-axis."""

    # throw out most diffuse functions
    alphas = [3.622534,6.274596,10.867601,18.823789]

    # ! assume proton is at the origin
    xo = yo = zo = 0

    # pick basis functions
    et_basis = []

    # add s-functions
    for alpha in alphas:
        et_basis.append({
            "l":0,
            "m":0,
            "n":0,
            "alpha":alpha,
            "xo":xo,
            "yo":yo,
            "zo":zo
        })
    # end for alpha

    # add p-functions
    for alpha in alphas:
        et_basis.append({
            "l":0,
            "m":0,
            "n":1,
            "alpha":alpha,
            "xo":xo,
            "yo":yo,
            "zo":zo
        })
    # end for alpha

    # add d-functions
    for alpha in alphas:
        et_basis.append({
            "l":2,
            "m":0,
            "n":0,
            "alpha":alpha,
            "xo":xo,
            "yo":yo,
            "zo":zo
        })
        
        et_basis.append({
            "l":0,
            "m":2,
            "n":0,
            "alpha":alpha,
            "xo":xo,
            "yo":yo,
            "zo":zo
        })
        
        et_basis.append({
            "l":0,
            "m":0,
            "n":2,
            "alpha":alpha,
            "xo":xo,
            "yo":yo,
            "zo":zo
        })
    # end for alpha

    return et_basis
# end def et20_basis
