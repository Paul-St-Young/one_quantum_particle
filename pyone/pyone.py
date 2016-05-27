#!/usr/bin/env python

# import python standard libraries
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

# import local libraries
import sys
sys.path.insert(0,"lib")
from basis import gto_dict,et20_basis
import grid


def total_energy(coeffs,basis_overlap_mat,basis_lapacian_mat
                 ,basis_potential_mat,mass=1836.):
    """ calculate kinetic and potential energies of the particle (default proton)
    only integral tables are needed. """

    nbasis = len(basis_overlap_mat)
    if len(coeffs) != nbasis:
        print "number of coeffcient %d != number of basis functions %d" % (len(coeffs),nbasis)
        return
    # end if
    
    # call wave function chi. It is a LC of basis functions
    chi_chi     = np.dot(coeffs, np.dot(basis_overlap_mat  ,coeffs) )
    chi_lap_chi = np.dot(coeffs, np.dot(basis_lapacian_mat ,coeffs) )
    chi_pot_chi = np.dot(coeffs, np.dot(basis_potential_mat,coeffs) )
    
    normalization    = chi_chi
    kinetic_energy   = -1./(2.*mass)*chi_lap_chi / normalization
    potential_energy = chi_pot_chi / normalization
    
    return kinetic_energy + potential_energy

# end def total_energy


if __name__ == "__main__":

    domain = grid.Potential_Grid()
    x,y,z,dat = domain.read_grid_ang("interp_grid/fat32.dat")
    print "contructing laplacian on %dx%dx%d grid" % (len(x),len(y),len(z))
    lap3d = domain.init_lap3d()

    et_basis = et20_basis()
    nbasis = len(et_basis)

    print "obtaining grid representations of basis functions"
    print " give me 2min ...  (for 32^3 grid)"
    # put basis functions down on a grid
    basis_set_grid = []
    for ibasis in range(nbasis):
        basis_set_grid.append( domain.grid_rep(
            lambda x,y,z:gto_dict(x,y,z,et_basis[ibasis])
        ) )
    # end for ibasis
    basis_set_grid = np.array( basis_set_grid )

    try:
        print "attemping to read integral tables"
        basis_overlap_mat   = np.loadtxt("int_table/et%d-overlap.dat"%nbasis)
        basis_lapacian_mat  = np.loadtxt("int_table/et%d-lapacian.dat"%nbasis)
        basis_potential_mat = np.loadtxt("int_table/et%d-potential.dat"%nbasis)
        print " read successful!"
    except:
        print "integral table read failed, constructing integral tables"
        # construct integral table for basis
        basis_overlap_mat  = np.zeros([nbasis,nbasis])
        basis_lapacian_mat = np.zeros([nbasis,nbasis])
        basis_potential_mat= np.zeros([nbasis,nbasis])
        for ibasis in range(nbasis):
            for jbasis in range(nbasis):
                
                if ibasis > jbasis:
                    # use symmetry to avoid redoing integral
                    basis_overlap_mat[ibasis,jbasis]   = basis_overlap_mat[jbasis,ibasis]
                    basis_lapacian_mat[ibasis,jbasis]  = basis_lapacian_mat[jbasis,ibasis]
                    basis_potential_mat[ibasis,jbasis] = basis_potential_mat[jbasis,ibasis]
                    continue
                # end if
                
                basis_overlap_mat[ibasis,jbasis] = sum(
                    basis_set_grid[ibasis]*basis_set_grid[jbasis] 
                ) * domain.measure
                
                basis_lapacian_mat[ibasis,jbasis] = sum(
                    basis_set_grid[ibasis]*np.dot( domain.lap3d,basis_set_grid[jbasis] )
                ) * domain.measure
                
                basis_potential_mat[ibasis,jbasis] = sum(
                    basis_set_grid[ibasis]*domain.dat*basis_set_grid[jbasis]
                ) * domain.measure
                
            # end for jbasis
        # end for ibasis

        # save integral tables
        np.savetxt("int_table/et%d-overlap.dat"%nbasis,basis_overlap_mat)
        np.savetxt("int_table/et%d-lapacian.dat"%nbasis,basis_lapacian_mat)
        np.savetxt("int_table/et%d-potential.dat"%nbasis,basis_potential_mat)
    # end try

    just_total = lambda guess: total_energy(guess
    ,basis_overlap_mat,basis_lapacian_mat,basis_potential_mat)

    guess = np.array([ 0.57619876,  0.62380963,  0.6189722 ,  0.33807449, -0.03288967,
            0.02785486,  0.21850042, -0.04154042, -0.03366639, -0.03364001,
           -0.05118077,  0.44831886,  0.44830779, -0.35853746,  0.32523168,
            0.32521758, -0.18786337,  0.01450228,  0.0145027 ,  0.04061571])

    print "optimizing coefficients, rerun with output as new guess if not converged"
    print 
    popt  = op.fmin(just_total,guess)
    
    # save optimized coefficients
    np.savetxt("popt.dat",popt)

    # save optimized wave function on the grid
    ground = np.dot( basis_set_grid.transpose(), popt )
    np.savetxt("ground.dat",ground)
    ground = ground.reshape(32,32,32)

    # extract bend and stretch 
    def gauss1d(x,A,xo,alpha):
        return A*np.exp(-alpha*(x-xo)**2.)
    # end def gauss1d
    xfit,xcov = op.curve_fit(gauss1d,domain.x,ground[:,16,16],p0=(5,0.1,4))
    zfit,zcov = op.curve_fit(gauss1d,domain.z,ground[16,16,:],p0=(5,0.1,4))
    print "bend    = %4.0f cm-1" % (xfit[-1]*237.)
    print "stretch = %4.0f cm-1" % (zfit[-1]*237.)

# end if
