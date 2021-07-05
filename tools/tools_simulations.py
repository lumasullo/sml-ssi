# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:29:23 2019

@author: Luciano Masullo, Lars Richter and Lucía López

Different tools used for MINFLUX exp analysis and simulations

"""
import numpy as np
import scipy.stats as stats

π = np.pi


def space_to_index(space, size_nm, px_nm):
    
    """
    
    NOTE: with the definiton of space_to_index and index_to_space, an array 
    with size_nm = 200 nm and px_nm will go from x = -100 nm to x = 99 nm and 
    y = -99 nm to y = 100 nm. This is arbitrary and ok as long as it is self 
    consistent, i.e. this functions are used to go from (x, y) to (i, j). All 
    space-coordinates array should also be consistent with this.
    
    """

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

def cov_ellipse(cov, q=None, nsig=None, **kwargs):
    
    """
    Plot of covariance ellipse
    
    Parameters
    ----------
    cov : (2, 2) array
        Covariance matrix.
    q : float, optional
        Confidence level, should be in (0, 1)
    nsig : int, optional
        Confidence level in unit of standard deviations. 
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.

    Returns
    -------
    width(w), height(h), rotation(theta in degrees):
         The lengths of two axises and the rotation angle in degree
    for the ellipse.
    """
    
    if q is not None:
        q = np.asarray(q)
    elif nsig is not None:
        q = 2 * stats.norm.cdf(nsig) - 1
    else:
        raise ValueError('Either `q` and `nsig` should be specified.')
    
    val, vec =  np.linalg.eig(cov)
    order = val.argsort()[::]
    val = val[order]
    vec = vec[order]
    w, h = 2 * np.sqrt(val[:, None])

    theta = np.degrees(np.arctan2(*vec[::, 0]))
    
    return w, h, theta


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
    
    P = 1
    d = 1.2*fwhm
    β = (4/π) * 2 * np.log(2) * (1/d**2)  # normalization to 1
    I = β * 2 * np.log(2) * (r**2/d**2) * np.exp(-4 * np.log(2) * (r**2/d**2))

    return I


def psf(central_zero, size, px, fwhm, psf_type):
    
    """ 2D extension for a 1D function in (r,θ) """

    x0 = central_zero[0]
    y0 = central_zero[1]
    
    x = np.arange(-size/2, size/2, px)
    y = -np.arange(-size/2, size/2, px) 
    
    # a minus sign is used in y to match cartesian convention
    # note that reversing the array y would imply a mismatch of 1 px with
    # respect to the index_to_space / space_to_index functions and thus
    # introduce an error
    
#    print('x', x[0], x[-1])
#    print('y', y[0], y[-1])

    
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

 
def ebp_centres(K, L, center, arr_type, phi=0):
    
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
        
        l = L / np.sqrt(2)
        
        pos_nm = np.zeros((K, 2))
        
        k = np.sqrt(K)
        
        # alternative definition of L
#        x_exc = np.arange(0, k) * L/k + (L/k)/2 - L/2
#        y_exc = np.arange(0, k)[::-1] * L/k + (L/k)/2 - L/2
        
        x_exc = np.arange(0, k) * l/(k-1) - l/2
        y_exc = np.arange(0, k)[::-1] * l/(k-1) - l/2

        # y_exc is reversed to match cartesian convention, (L/K)/2 is added to place
        
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


def sim_exp(psf, r0_nm, N, sbr, size_nm, px_nm, localizations=1, DEBUG=False):
    
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
    
    r0 = space_to_index(r0_nm, size_nm, px_nm)

    K = np.shape(psf)[0] # number of expositions
            
    # normalization term, sum of all psf intensities
    norm_psf = np.sum(psf, axis=0)[r0[0], r0[1]]
    
    # p-parameter for each beam
    
    p = np.zeros(K)
    
    for i in range(K):
        
        p[i] = sbr/(sbr+1) * psf[i, r0[0], r0[1]]/norm_psf + (1/(sbr+1)) * (1/K)
            
    if localizations != 1:    
        
        n_array = np.zeros((localizations, K))
        
        for i in range(localizations):
                    
            n_array[i, :] = np.random.multinomial(N, p)
        
    else:
        
        n_array = np.random.multinomial(N, p)
        
    return n_array


def crb(K, psf, sbr, px_nm, size_nm, N, prior='rough loc', s=50):
    
    """
    
    Cramer-Rao Bound for a given SSI-SML experiment 
    
    Input
    ----------
    K : int, number of excitation beams
    
    psf: (K, size, size) array, experimental or simulated excitation beams
    
    SBR : float, signal to background ratio
    
    px_nm : pixel of the grid in nm
    
    size_nm : size of the grid in nm
    
    N : total number of photons
    
    prior: type of a priori information, set to None if no prior is included
    'rough loc' is a model for a previous rough localization (gaussian-like)
    
    s: σ parameter of the gaussian information prior, s=50 nm as default
    rough localization, s >> L would be equivalent to no prior, s << L might
    not be a realistic assumption
    
    Output: σ_CRB, Σ_CRB, Fr, sbr_rel
    
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
    
    sbr_rel = np.sum(psf, axis=0)/np.sum(psf, axis=0)[int(size/2), int(size/2)]
    
    # SBR is computed as SBR(x, y), i.e. proportional to rel_sbr(x, y) instead 
    # of a constant
    
    sbr = sbr_rel * sbr
    
    # normalization term, sum of all psf intensities
    norm_psf = np.sum(psf, axis = 0)
    
    for i in range(K):
        
        p[i, :, :] = sbr/(sbr+1) * psf[i,:,:]/norm_psf + (1/(sbr+1)) * (1/K)
        
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
                                                        
                Fr[:, :, i, j] = N * np.sum(Fr_aux[:, :, :, i, j], axis=0) + np.diag([1/s**2, 1/s**2])
                            
            elif prior is None:
                          
                Fr[:, :, i, j] = N * np.sum(Fr_aux[:, :, :, i, j], axis=0)
                
            Σ_CRB[:, :, i, j] = np.linalg.inv(Fr[:, :, i, j])
            σ_CRB[i, j] = np.sqrt((1/d) * np.trace(Σ_CRB[:, :, i, j]))
            
    return σ_CRB, Σ_CRB, Fr, sbr_rel

   
def mle(n, psf, sbr, px_nm=1, prior=None, s=None, localizations=1, 
        DEBUG=False):
    
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
    
    #TO DO: fix for px_nm != 1
       
    # number of beams in EBP
    K = np.shape(psf)[0]
    
    # FOV size
    size = np.shape(psf)[1] 
    
    # normalization term, sum of all psf intensities
    norm_psf = np.sum(psf, axis = 0)

    # sbr(x, y) is computed as SBR(x, y), i.e. proportional to rel_sbr(x, y) 
    # instead of a constant
    sbr_rel = np.sum(psf, axis=0)/np.sum(psf, axis=0)[int(size/2), int(size/2)]
    sbr = sbr_rel * sbr
    
    # p-parameter vector 
    p = np.zeros((K, size, size))
    
    for i in np.arange(K):        
        p[i,:,:] = sbr/(sbr + 1) * psf[i,:,:]/norm_psf + (1/(sbr + 1)) * (1/K)
    
    if localizations != 1:
        
        mle_index = np.zeros((localizations, 2))
        mle_nm = np.zeros((localizations, 2))
        likelihood = np.zeros((localizations, size, size))
        
        for i in range(localizations):
            
            print('localization ' + str(i))
            
            # log-likelihood function
            l_aux = np.zeros((K, size, size))
            for k in range(K):
                l_aux[k, :, :] = n[i, k] * np.log(p[k, : , :])
                
            likelihood[i, :, :] = np.sum(l_aux, axis=0)
                        
            if prior == 'rough loc':
                                            
                x = np.arange(-size/2, size/2) * px_nm
                y = np.arange(-size/2, size/2) * px_nm
                
                xx, yy = np.meshgrid(x, y)
                                
                # log of 2D gaussian function
                prior_info = -(xx**2)/(2*s**2)-(yy**2)/(2*s**2)
                
                likelihood[i, :, :] = likelihood[i, :, :] + prior_info
                    
            else:
            
                pass
            
            # maximum likelihood estimator for the position   
        
            mle_index[i, :] = np.unravel_index(np.argmax(likelihood[i, :, :], 
                                               axis=None), 
                                               likelihood[i, :, :].shape)
        
            mle_nm[i, :] = index_to_space(mle_index[i, :], size * px_nm, px_nm)
        
    else:
    
        # log-likelihood function
        l_aux = np.zeros((K, size, size))
        for i in range(K):
            l_aux[i, :, :] = n[i] * np.log(p[i, : , :])
            
        likelihood = np.sum(l_aux, axis=0)
                        
        if prior == 'rough loc':
                                        
            x = np.arange(-size/2, size/2) * px_nm
            y = np.arange(-size/2, size/2) * px_nm
            
            xx, yy = np.meshgrid(x, y)
                            
            # log of 2D gaussian function
            prior_info = -(xx**2)/(2*s**2)-(yy**2)/(2*s**2)
            
            likelihood = likelihood + prior_info
                    
        else:
            
            pass
                             
        # maximum likelihood estimator for the position    
        mle_index = np.unravel_index(np.argmax(likelihood, axis=None), 
                                     likelihood.shape)
        
        mle_nm = index_to_space(mle_index, size * px_nm, px_nm)
    
    
    if DEBUG:
        
        return mle_nm, mle_index, likelihood
    
    else:   
    
        return mle_nm
    
    
def crb_camera(K, px_size_nm, sigma_psf, sbr, N, range_nm, dx_nm):
 
    """
    Cramer-Rao Bound for camera-based localization
    
    Input
    ----------
        K:          total number of physical pixels in grid,
                    e.g. K=81 if 9x9 grid
        px_size_nm: pixel size of camera in nm
        sigma_psf:  standard deviation of emission PSF in nm
        sbr:        signal to background ratio
        N :         total number of photons
        range_nm:   grid range in emitter position space
        dx_nm:      grid resolution along x-axis in emitter position space
    
    Returns
    -------
        crb:        Scalar CRB for the given set of input parameters
    
    """
    
    from scipy.special import erf
    
    σ_psf = sigma_psf
    dy_nm = dx_nm
    size = int(range_nm/dx_nm)
    
    if np.sqrt(K).is_integer(): 
        k = int(np.sqrt(K)) # k^2 is equivalent to the K expositions in MINFLUX CRB
    else:
        raise ValueError("K should be a perfect square number (i.e. 81, 64, etc)")
    
    p, λ, dpdx, dpdy, A, B, C, D = (np.zeros((size, size, k, k)) for i in range(8))
    
    x = np.arange(-range_nm/2, range_nm/2, dx_nm)
    y = np.arange(-range_nm/2, range_nm/2, dy_nm)
    
    px_range_nm = k * px_size_nm
    px = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    py = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    
    # -y for cartesian coordinates
    [Mx, My] = np.meshgrid(x, -y)  # emitter position matrices, dim: size x size (i.e. 100 x 100)
    [Mpx, Mpy] = np.meshgrid(px, -py) # pixel coord matrices dim: k x k (i.e. 9 x 9)
        
    for i in range(k):
        for j in range(k):
            
            # calculates terms
            a = erf((Mpx[i, j] + px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            b = erf((Mpx[i, j] - px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            c = erf((Mpy[i, j] + px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            d = erf((Mpy[i, j] - px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            
            p_0 = (1/4) * (a-b) * (c-d)
            
            # p array for a given pixel (i, j) and for every position (x, y)
            p[:, :, i, j] = sbr/(sbr+1) * p_0 + + (1/(sbr+1)) * (1/K)
            
            # gradients of p in each (x,y), careful with (i, j) -> (x, y)
            dpdy[:, :, i, j] = np.gradient(p[:, :, i, j], -dy_nm, axis=0)
            dpdx[:, :, i, j] = np.gradient(p[:, :, i, j], dx_nm, axis=1)
       
            # terms needed to compute CRB for each (x,y) in emitter space
            A[:, :, i, j] = (1/p[:, :, i, j]) * dpdx[:, :, i, j]**2
            B[:, :, i, j] = (1/p[:, :, i, j]) * dpdy[:, :, i, j]**2
            C[:, :, i, j] = (1/p[:, :, i, j]) *(dpdx[:, :, i, j] * dpdy[:, :, i, j])
            D[:, :, i, j] = (1/p[:, :, i, j]) * (dpdx[:, :, i, j]**2 + dpdy[:, :, i, j]**2)
                
    # sigma CRB analytical expression for the scalar form    
    E = np.sum(D, axis=(2,3)) 
    F = (np.sum(A, axis=(2,3)) * np.sum(B, axis=(2,3))) - np.sum(C, axis=(2,3))**2
    
    σ_CRB = np.sqrt(1/(2*N))*np.sqrt(E/F)
    
    return σ_CRB
    

def sim_exp_camera(r0_nm, K, px_size_nm, sigma_psf, sbr, N, range_nm, dx_nm):
    
    from scipy.special import erf
    
    σ_psf = sigma_psf
    dy_nm = dx_nm
    size = int(range_nm/dx_nm)
    
    if np.sqrt(K).is_integer(): 
        k = int(np.sqrt(K)) # k^2 is equivalent to the K expositions in MINFLUX CRB
    else:
        raise ValueError("K should be a perfect square number (i.e. 81, 64, etc)")
    
    p, λ, dpdx, dpdy, A, B, C, D = (np.zeros((size, size, k, k)) for i in range(8))
    
    x = np.arange(-range_nm/2, range_nm/2, dx_nm)
    y = np.arange(-range_nm/2, range_nm/2, dy_nm)
    
    px_range_nm = k * px_size_nm
    px = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    py = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    
    # -y for cartesian coordinates
    [Mx, My] = np.meshgrid(x, -y)  # emitter position matrices, dim: size x size (i.e. 100 x 100)
    [Mpx, Mpy] = np.meshgrid(px, -py) # pixel coord matrices dim: k x k (i.e. 9 x 9)
            
    for i in range(k):
        for j in range(k):
            
            # calculates terms
            a = erf((Mpx[i, j] + px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            b = erf((Mpx[i, j] - px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            c = erf((Mpy[i, j] + px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            d = erf((Mpy[i, j] - px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            
            p_0 = (1/4) * (a-b) * (c-d)
            
            # p array for a given pixel (i, j) and for every position (x, y)
            p[:, :, i, j] = sbr/(sbr+1) * p_0 + + (1/(sbr+1)) * (1/K)
    
    r0 = space_to_index(r0_nm, range_nm, dx_nm)        
    p_r0 = p[r0[0], r0[1], :, :].ravel()
        
    n_array = np.random.multinomial(N, p_r0)
    
    return n_array    
    
    
def mle_camera(n, K, px_size_nm, sigma_psf, sbr, N, range_nm, dx_nm):
    
    from scipy.special import erf
    
    σ_psf = sigma_psf
    dy_nm = dx_nm
    size = int(range_nm/dx_nm)
    
    if np.sqrt(K).is_integer(): 
        k = int(np.sqrt(K)) # k^2 is equivalent to the K expositions in MINFLUX CRB
    else:
        raise ValueError("K should be a perfect square number (i.e. 81, 64, etc)")
    
    p, λ, dpdx, dpdy, A, B, C, D = (np.zeros((size, size, k, k)) for i in range(8))
    
    x = np.arange(-range_nm/2, range_nm/2, dx_nm)
    y = np.arange(-range_nm/2, range_nm/2, dy_nm)
    
    px_range_nm = k * px_size_nm
    px = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    py = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
    
    # -y for cartesian coordinates
    [Mx, My] = np.meshgrid(x, -y)  # emitter position matrices, dim: size x size (i.e. 100 x 100)
    [Mpx, Mpy] = np.meshgrid(px, -py) # pixel coord matrices dim: k x k (i.e. 9 x 9)
        
    for i in range(k):
        for j in range(k):
            
            # calculates terms
            a = erf((Mpx[i, j] + px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            b = erf((Mpx[i, j] - px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
            c = erf((Mpy[i, j] + px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            d = erf((Mpy[i, j] - px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
            
            p_0 = (1/4) * (a-b) * (c-d)
            
            # p array for a given pixel (i, j) and for every position (x, y)
            p[:, :, i, j] = sbr/(sbr+1) * p_0 + + (1/(sbr+1)) * (1/K)
            
            
    p = p.reshape(size, size, K)        
    # log-likelihood function
    l_aux = np.zeros((K, size, size))
    for i in np.arange(K):
        l_aux[i, :, :] = n[i] * np.log(p[: , :, i])
        
    likelihood = np.sum(l_aux, axis = 0)
    
    # maximum likelihood estimator for the position   
    
    mle_index = np.unravel_index(np.argmax(likelihood, axis=None), 
                                 likelihood.shape)
    
    mle_nm = index_to_space(mle_index, size, dx_nm)

    return mle_nm


#def crb_camera(K, px_size_nm, sigma_psf, SBR, N, range_nm, dx_nm):
# 
#    """
#    Cramer-Rao Bound for camera-based localization
#    
#    Input
#    ----------
#        px_num:     total number of physical pixels in grid,
#                    e.g. px_num=81 if 9x9 grid
#        px_size_nm: pixel size of camera in nm
#        sigma_psf:  standard deviation of emission PSF in nm
#        SBR:        signal to background ratio
#        N :         total number of photons
#        range_nm:   grid range in emitter position space
#        dx_nm:      grid resolution along x-axis in emitter position space
#        inf:        boolean, if true, set SBR to infinity, default: false
#    
#    Returns
#    -------
#        crb:        Cramer-Rao bound for given set of input parameter
#    
#    """
#    
#    from scipy.special import erf
#    
#    σ_psf = sigma_psf
#    dy_nm = dx_nm
#    size = int(range_nm/dx_nm)
#    
#    if np.sqrt(K).is_integer(): 
#        k = int(np.sqrt(K)) # k^2 is equivalent to the K expositions in MINFLUX CRB
#    else:
#        raise ValueError("K should be a perfect square number (i.e. 81, 64, etc)")
#    
#    p, λ, dpdx, dpdy, A, B, C, D = (np.zeros((size, size, k, k)) for i in range(8))
#    
#    x = np.arange(-range_nm/2, range_nm/2, dx_nm)
#    y = np.arange(-range_nm/2, range_nm/2, dx_nm)
#    
#    px_range_nm = k * px_size_nm
#    px = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
#    py = np.arange(-px_range_nm/2, px_range_nm/2, px_size_nm) + px_size_nm/2
#    
#    # -y for cartesian coordinates
#    [Mx, My] = np.meshgrid(x, -y)  # emitter position matrices, dim: size x size (i.e. 100 x 100)
#    [Mpx, Mpy] = np.meshgrid(px, -py) # pixel coord matrices dim: k x k (i.e. 9 x 9)
#        
#    for i in range(k):
#        for j in range(k):
#            
#            # calculates terms
#            a = erf((Mpx[i, j] + px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
#            b = erf((Mpx[i, j] - px_size_nm/2 - Mx)/(np.sqrt(2)*σ_psf))
#            c = erf((Mpy[i, j] + px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
#            d = erf((Mpy[i, j] - px_size_nm/2 - My)/(np.sqrt(2)*σ_psf))
#            
#            p_0 = (a-b) * (c-d)
#            
#            # prob array for a given pixel (i, j) and for every position (x, y)
#            p[:, :, i, j] = (1/(K + SBR)) + (1/4) * (SBR/(K + SBR)) * p_0
#            
#            # gradient of ps in each (x,y), careful with (i, j) -> (x, y)
#            dpdy[:, :, i, j] = np.gradient(p[:, :, i, j], -dy_nm, axis=0)
#            dpdx[:, :, i, j] = np.gradient(p[:, :, i, j], dx_nm, axis=1)
#       
#            # terms needed to compute CR bound for each (x,y) in emitter space
#            A[:, :, i, j] = (1/p[:, :, i, j]) * dpdx[:, :, i, j]**2
#            B[:, :, i, j] = (1/p[:, :, i, j]) * dpdy[:, :, i, j]**2
#            C[:, :, i, j] = (1/p[:, :, i, j]) *(dpdx[:, :, i, j] * dpdy[:, :, i, j])
#            D[:, :, i, j] = (1/p[:, :, i, j]) * (dpdx[:, :, i, j]**2 + dpdy[:, :, i, j]**2)
#                
#    # sigma Cramer-Rao numerator and denominator    
#    E = np.sum(D, axis=(2,3)) 
#    F = (np.sum(A, axis=(2,3)) * np.sum(B, axis=(2,3))) - np.sum(C, axis=(2,3))**2
#    
#    σ_CRB = np.sqrt(1/(2*N))*np.sqrt(E/F)
#    
#    return σ_CRB
    
