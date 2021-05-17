# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:29:23 2019

@author: Luciano Masullo, Lars Richter and Lucía López

Different tools used for MINFLUX exp analysis and simulations

"""
import numpy as np

π = np.pi


def space_to_index(space, size_nm, px_nm):

    # size and px have to be in nm
    index = np.zeros(2)
    index[0] = (size_nm/2 - space[1])/px_nm
    index[1] = (space[0]+ size_nm/2)/px_nm 

    return np.array(index, dtype=np.int)


def index_to_space(index, size_nm, px_nm):

    space = np.zeros(2)
    space[0] = index[1]*px_nm - size_nm/2
    space[1] = size_nm/2 - index[0]*px_nm

    return np.array(space)


def create_circular_mask(h, w, center=None, radius=None):
    
    """    
    Auxiliar function to create a circular mask for an array
    """    

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def gaussian(r, fwhm):
    
    """ 2D gaussian beam intensity """
    
    β = 1
    I = β * np.exp(-4 * np.log(2) * (r**2/fwhm**2))
    
    return I


def doughnut(r, fwhm):
    
    """ 2D donut """

    d = 1.2*fwhm
#    β = (4/π) * 2 * np.log(2) * (P/d**2)  # normalization to 1
    β = 2*np.e # TODO: check normalization
    I = β * 2 * np.log(2) * (r**2/d**2) * np.exp(-4 * np.log(2) * (r**2/d**2))

    return I


def psf(central_zero, size, px, fwhm, psf_type):
    
    """ 2D extension for a 1D function in (r,θ) """

    x0 = central_zero[0]
    y0 = central_zero[1]

    x = np.arange(-size/2, size/2, px)
    y = np.arange(-size/2, size/2, px)[::-1] # to match cartesian convention

    [Mx, My] = np.meshgrid(x, y)
    
    if psf_type == 'doughnut':
        Mro = np.sqrt((Mx-x0)**2 + (My-y0)**2)
        result = doughnut(Mro, fwhm)
    
    elif psf_type == 'gaussian':
        Mro = np.sqrt((Mx-x0)**2 + (My-y0)**2)
        result = gaussian(Mro, fwhm)
        
    else:
        
        print('Please select a valid psf. Valid types: gaussian, doughnut')
    
    return result

 
def ebp_centres(K, L, center, phi=0, arr_type='orbit'):
    
    """   
    Calculates doughnut centre positions of the excitation pattern.

    Input
    ----------
        K: number of distinct doughnuts, usually K=4
        L: EBP diameter
        center: boolean, true/false: EBP with or without centre
        arr_type: orbit or raster scan type of EBP
    
    Returns
    -------
        pos_nm: EBP beam centre positions
       
    """
    
    if arr_type == 'raster scan':
        
        pos_nm = np.zeros((K, 2))
        
        k = np.sqrt(K)
        x_exc = np.arange(0, k) * L/k + (L/k)/2 - L/2
        y_exc = np.arange(0, k)[::-1] * L/k + (L/k)/2 - L/2

        # y_exc is reversed to match cartesian convention, (L/K)/2 is added to place
        # the excitation position in the center of the step (px in acquired image)
        # i.e. [  10.  30.  50.  70.  90. 110.] for a L = 120 nm configuration
        
        xx_exc, yy_exc = np.meshgrid(x_exc, y_exc)

        pos_nm[:, 0] = xx_exc.ravel()
        pos_nm[:, 1] = yy_exc.ravel()
    
    elif arr_type == 'orbit':
    
        pos_nm = np.zeros((K, 2)) 
        θ = np.zeros(K)
        L = np.ones(K)*L
    #    θrandom = [π/10, π/12.3, -π/18, -π/9]
    
        if center:
            # beams centers (for regular polygons incribed in a circle of diameter L)
            pos_nm[0, :] = 0. # central doughnut at (0, 0)
            Kθ = K - 1
            if (Kθ+1) % 2 == 0:  # if K is odd
                for k in np.arange(1, Kθ+1):
                    θ[k] = π * 2*k/Kθ + phi 
                    pos_nm[k, 0] = (L[k]/2) * np.cos(θ[k])
                    pos_nm[k, 1] = (L[k]/2) * np.sin(θ[k])
            else:       # if K is even
                for k in np.arange(1, Kθ+1):
                    θ[k] = π * (2*k+1)/Kθ + phi 
                    pos_nm[k, 0] = (L[k]/2) * np.cos(θ[k])
                    pos_nm[k, 1] = (L[k]/2) * np.sin(θ[k])
        else:
            if (K+1) % 2 == 0:  # if K is odd
                for k in np.arange(K):
                    θ[k] = π * 2*k/K + phi 
                    pos_nm[k, 0] = (L[k]/2) * np.cos(θ[k])
                    pos_nm[k, 1] = (L[k]/2) * np.sin(θ[k])
            else:       # if K is even
                for k in np.arange(K):
                    θ[k] = π * (2*k+1)/K + phi 
                    pos_nm[k, 0] = (L[k]/2) * np.cos(θ[k])
                    pos_nm[k, 1] = (L[k]/2) * np.sin(θ[k])
                    
    else:
        
        print('Please choose a valid type of EBP')
                        
    return pos_nm  


def sim_exp(psf, r0, N, SBR, DEBUG=False):
    
    """
    
    This function generates a simulated measurement of an ssi-sml experiment.
    
    Input:
    
    - psf: np.array (float)
    
        Stack of K PSFs of the minflux experiment
    
    - r0: np.array (float)
    
        2D position of the emitter (x, y)
        
    - N: int
        
        Number of total detected photons
        
    - SBR: int
    
        Signal to background ratio
        
    Output:
        
    - n_array: np.array (int)
    
        Array of measured photons at each exposition
    
    """

    K = np.shape(psf)[0] # number of expositions
        
    λ = np.zeros(K)
    for i in np.arange(K):
        λ[i] = psf[i, r0[0], r0[1]] 
    
    # backgroung given a SBR level
    λb = np.sum(λ)/(K*SBR)
    normλ = (np.sum(λ) + K*λb)
    
    # probability for each beam
    p = np.zeros(K)
    for i in np.arange(K):
        p[i] = (λ[i] + λb)/ normλ
                    
    n_array = np.random.multinomial(N, p)
    
    return n_array
    

def crb(K, psf, SBR, px_nm, size_nm, N, prior='rough loc'):
    
    """
    
    Cramer-Rao Bound for a given SSI-SML experiment 
    
    Input
    ----------
    K : int, number of excitation beams
    psf_array : (K, size, size) array, experimental or simulated PSF 
    SBR : float, signal to background ratio
    px_nm : pixel of the grid in nm
    size_nm : size of the grid in nm
    N : total number of photons
    
    Method: calculates the Σ_CRB from the Fisher information matrix in 
    emitter position space (Fr), from there it calculates Σ_CRB and σ_CRB
    (S11-13, 10.1126/science.aak9913)
    
    Output: Fr, Σ_CRB, σ_CRB
    
    """
    
    # size of the σ_CRB matrix in px and dimension d=2
    size = int(size_nm/px_nm)
    d = 2
    
    # size of the (x,y) grid
    dx = px_nm
    dy = px_nm
        
    # initialize different arrays needed to compute σ_CRB, Σ_CRB and Fr
    
    σ_CRB = np.zeros((size, size))
    p, λ, dpdx, dpdy = (np.zeros((K, size, size)) for i in range(4))
    Fr, Σ_CRB = (np.zeros((d, d, size, size)) for i in range(2))
    
    Fr_aux = np.zeros((K, d, d, size, size))
    
    # normalization of PSF to Ns = N*(SBR/(SBR+1))

    for i in range(K):
        
        λ[i, :, :] = N*(SBR/(SBR+1)) * (psf[i, :, :]/np.sum(psf, axis=0))
        
    # λb using the approximation in Balzarotti et al, (S29)
      
    λb = np.sum(λ[:, int(size/2), int(size/2)])/(K*SBR)
        
    for i in range(K):
        
        # p-parameter arrays
    
        p[i, :, :] = (λ[i, :, :] + λb)/(K*λb + np.sum(λ, axis=0))

        # partial derivatives in x and y (minus sign for cartesian convention)

        dpdy[i, :, :], dpdx[i, :, :] = np.gradient(p[i, :, :], -dy, dx)
        
    # compute relevant information for every (i, j) position
    
    # TODO: vectorize this part of the code
            
    for i in range(size):
        for j in range(size):
            
            for k in range(K):
        
                A = np.array([[dpdx[k, i, j]**2, 
                               dpdx[k, i, j]*dpdy[k, i, j]],
                              [dpdx[k, i, j]*dpdy[k, i, j], 
                               dpdy[k, i, j]**2]])
    
                Fr_aux[k, :, :, i, j] = (1/p[k, i, j]) * A
                
            if prior == 'rough loc':
                            
                s = 50 # in nm, sigma of the prior rough gaussian localization
                            
                Fr[:, :, i, j] = N * np.sum(Fr_aux[:, :, :, i, j], axis=0) + np.diag([1/s**2, 1/s**2])
                            
            elif prior is None:
                          
                Fr[:, :, i, j] = N * np.sum(Fr_aux[:, :, :, i, j], axis=0)
                
            Σ_CRB[:, :, i, j] = np.linalg.inv(Fr[:, :, i, j])
            σ_CRB[i, j] = np.sqrt((1/d) * np.trace(Σ_CRB[:, :, i, j]))
            
    return σ_CRB, Σ_CRB, Fr

     
def mle(n, PSF, SBR, px_nm=1, prior=None, s=None, DEBUG=False):
    
    """    
    Position estimator (using MLE)
    
    Inputs
    ----------
    n : acquired photon collection (K)
    PSF : array with EBP (K x size x size)
    SBR : estimated (exp) Signal to Background Ratio
    prior: prior information on the emitter position
    s: parameter of the prior (sigma of the rough gaussian loc or radius of the 
                               search area)

    Returns
    -------
    mle_index : position estimator in index coordinates (MLE)
    likelihood : Likelihood function
    
    Parameters 
    ----------
    px_nm : grid px in nm
        
    """
       
    # number of beams in EBP
    K = np.shape(PSF)[0]
    
    # FOV size
    size = np.shape(PSF)[1] 
    
    normPSF = np.sum(PSF, axis = 0)
    
    # probabilitiy vector 
    p = np.zeros((K, size, size))

    for i in np.arange(K):        
        p[i,:,:] = (SBR/(SBR + 1)) * PSF[i,:,:]/normPSF + (1/(SBR + 1)) * (1/K)
        
    # log-likelihood function
    l_aux = np.zeros((K, size, size))
    for i in np.arange(K):
        l_aux[i, :, :] = n[i] * np.log(p[i, : , :])
        
    likelihood = np.sum(l_aux, axis = 0)
        
    if prior == 'r<s':
                
        x = np.arange(-size/2, size/2)
        y = np.arange(-size/2, size/2)
        
        Mx, My = np.meshgrid(x, y)
        Mr = np.sqrt(Mx**2 + My**2)
        
        likelihood[Mr>s/2] = -np.inf
        
    elif prior == 'rough loc':
                        
        x = np.linspace(-size/2, size/2, size)
        y = np.linspace(-size/2, size/2, size)
        
        xx, yy = np.meshgrid(x, y)
                
        G = np.exp(-(xx**2)/(s**2)-(yy**2)/(s**2))
        
        likelihood = likelihood + np.log(G)
        
    else:
        
        pass
                         
    # maximum likelihood estimator for the position   
    
    mle_index = np.unravel_index(np.argmax(likelihood, axis=None), 
                                 likelihood.shape)
    
    mle_nm = index_to_space(mle_index, size, px_nm)
    
    if DEBUG:
        
        return mle_nm, mle_index, likelihood
    
    else:   
    
        return mle_nm
