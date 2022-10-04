# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:07:30 2019

@author: simkon
"""

import numpy as np

def matCalc_concrete(eps0, deps, time, dtime, predef, dpredef, matProps, stress, stiff, stateVars):
    
    E_0, ft, gamma, fc, eps_c0, alpha_critical, sig_u, sig_L, beta_E = matProps[:9]
    E_creep = matProps[9:14]
    lambda_creep_array = matProps[14:19]
    beta_Ecreep = np.array([[-0.01914656,  0.07135232,  0.26914982,  0.227365  ,  0.11070075],
                            [-0.05485107,  0.11590657,  0.30621064,  0.26957512,  0.27000593],
                            [-0.05014045,  0.1123435 ,  0.24394659,  0.23264344,  0.18283563],
                            [-0.0339299 ,  0.09002191,  0.28048416,  0.23998737,  0.18048002],
                            [-0.0555953 ,  0.11949636,  0.24336659,  0.24786201,  0.15465828]])
    stress, stiff, stateVars = np.copy(stress), np.copy(stiff), np.copy(stateVars)
    alpha, eps_asr, alphac = stateVars[:3]
    alphac = max(alphac,0.000001)
    zeta, xi = stateVars[4:9], stateVars[9:14]

    N_creep = len(E_creep)

    
    
    #--STRAIN CALCULATIONS--
    deps_asr =   depsAsrCalc(stress, dpredef, sig_u, sig_L) 
    eps_asr = eps_asr + deps_asr
    
    eps_creep = epsCreepCalc(lambda_creep_array, zeta, xi, dtime)
    
    eps = eps0 + deps
    eps_i = eps - eps_asr - eps_creep
    
    
    sig = stiff*eps_i
    Ds = Ds_calc(E_0, ft, fc, eps_c0, 0., beta_E, eps_asr, alpha, alphac, sig)
    sig = Ds * eps_i
    
    
    
    sigcr,dsigcr_dalpha = sigcr_calc(E_0, ft, 0., beta_E, eps_asr, alpha, alpha_critical)
    fcr = sig - sigcr
    
    sigc,dsigc_dalphac = sigc_calc(fc, E_0, eps_c0, alphac)
    
    Fc = -sig - sigc
    if fcr > 0.:
        N_CONVERGENCE = 0
        for itr in range(20):
            dCcr_dalpha = (sigcr - alpha*dsigcr_dalpha)/(sigcr**2.)
            dfcr_dalpha = -Ds*dCcr_dalpha*sig - dsigcr_dalpha
            dalpha = -fcr/dfcr_dalpha
            alpha = alpha + dalpha
            alpha = max(alpha, stateVars[0])
            Ds = Ds_calc(E_0, ft, fc, eps_c0, 0., beta_E, eps_asr, alpha, alphac, sig)
            sig = Ds * eps_i
            sigcr,dsigcr_dalpha = sigcr_calc(E_0, ft, 0., beta_E, eps_asr, alpha, alpha_critical)
			
            fcr = sig - sigcr
            if (abs(fcr) <= 0.000001):
                N_CONVERGENCE = 1
                break
        if (N_CONVERGENCE != 1):
            print('DID NOT CONVERGE ON MATERIAL LEVEL')
    elif Fc > 0.:
        N_CONVERGENCE = 0
        for itr in range(20):
            dCdalphac = (sigc - alphac*dsigc_dalphac)/(sigc**2.)
            dFcdalphac = Ds*dCdalphac*sig - dsigc_dalphac
            dalphac = -Fc/dFcdalphac
            alphac = alphac + dalphac
            alphac = max(alphac, stateVars[2])
            Ds = Ds_calc(E_0, ft, fc, eps_c0, 0., beta_E, eps_asr, alpha, alphac, sig)
            sig = Ds * eps_i
            sigc,dsigc_dalphac = sigc_calc(fc, E_0, eps_c0, alphac)
			
            Fc = -sig - sigc
            if (abs(Fc) <= 0.000001):
                N_CONVERGENCE = 1
                break
        if (N_CONVERGENCE != 1):
            print('DID NOT CONVERGE ON MATERIAL LEVEL')
        

    stiff =  Ds
    dsig = sig - stress
    stress  = sig
    stateVars[:3] = alpha, eps_asr, alphac
    
    time_mid = time + dtime/2.
    for k in range(N_creep):
        Ccreep_k = creepCompliance_calc(lambda_creep_array, E_creep[k], beta_Ecreep[k], time_mid)
        zeta[k] = zeta[k] + Ccreep_k *dsig
        xi[k]   = np.exp(-dtime/lambda_creep_array[k]) * xi[k] + lambda_creep_array[k]*Ccreep_k/dtime*(1-np.exp(-dtime/lambda_creep_array[k]))*dsig
    stateVars[4:9], stateVars[9:14] = zeta, xi

    return stress, stiff, stateVars

def matCalc_steel(eps0, deps, time, dtime, predef, dpredef, matProps,sig, stiff,stateVars):
    

    E_0, f_y, sh = matProps 
    sig, stiff, stateVars = np.copy(sig), np.copy(stiff), np.copy(stateVars)
    eps_pl, kappa = np.copy(stateVars)
    
    eps = eps0 + deps
    sig = E_0 * (eps-eps_pl)
    
    sig_y = f_y + sh * E_0 *kappa
    f_p = sig - sig_y
    if f_p > 0.:
        #plasticity
        #print('Yielding')
        n_convergence = False
        dkappa = 0.
        for i in range(10):
            ddkappa = 1/(E_0+sh*E_0)*f_p
            dkappa = dkappa + ddkappa
            kappa = kappa + ddkappa
            sig_y = f_y + sh * E_0 *kappa
            eps_pl = eps_pl + ddkappa * np.sign(sig)
            sig = E_0 *(eps-eps_pl)
            f_p = sig - sig_y
            if f_p < 10**-6:
                n_convergence = True
                stiff = sh*E_0/(1+sh)
                break
        if n_convergence == False:
            print('No convergence on material level')
    else:
        stiff = E_0
    
    stateVars = np.array([eps_pl, kappa])
    
        
    return sig, stiff, stateVars


def depsAsrCalc(sig, dpredef, sig_u, sig_L):
    
    #-----------ASR CALCULATIONS---------------
    
    deps_asr_free = dpredef

    if sig < sig_u:
        W = 0.
    elif sig < sig_L:
        W = 1. - np.log10(sig/sig_L)/np.log10(sig_u/sig_L)
    else:
        W = 1.
    if W < 0.:
        print('ERROR: NEGATIVE W')

    deps_asr  = deps_asr_free * W
    
    return deps_asr


def epsCreepCalc(lambda_array, zeta, xi, dtime):
    N_creep = len(lambda_array)
    eps_creep = 0.
    for k in range(N_creep):
        eps_creep = eps_creep + zeta[k] - np.exp(-dtime/lambda_array[k])*xi[k]
    
    return eps_creep



def creepCompliance_calc(lambda_array, Ecreep_k, beta_Ecreep_k, tau):
    
    N_creep = len(lambda_array)
      
    Ccreep_k = 1.0
    for j in range(N_creep):
        Ccreep_k = Ccreep_k - beta_Ecreep_k[j] * (1. - np.exp(-tau/lambda_array[j]) )

    Ccreep_k = Ccreep_k * 1./Ecreep_k
    return  Ccreep_k
	     
    
def sigcr_calc(E_0, ft, gamma, beta_E, eps_asr, alpha, alpha_critical):

    if alpha ==0.:#<= alpha_critical:
        sigcr = ft#*np.exp(-gamma*alpha/ft*E_0)
        dsigcr_dalpha = 0.#-ft*np.exp(-gamma*alpha/ft*E_0)*gamma/ft*E_0
    else:
        sigcr = 0.01*ft
        dsigcr_dalpha = 0.0
        #print('fully open crack')
		
    return sigcr, dsigcr_dalpha
	    
def sigc_calc(fc,Ec,eps_c0,alphac):
    m = eps_c0*Ec/fc
    n = 1./(1.-1./m)

    if np.abs(alphac) <= eps_c0:
        k = 1.
    else:
        k = 0.67 + fc/62
    
    if alphac < 0.0035:
        sigc = (n * fc)/(eps_c0*(n-1+(abs(alphac)/eps_c0)**(n*k))) * alphac
        dalphac = 0.00001
        dsigc_dalphac = ((n * fc)/(eps_c0*(n-1+(abs(alphac+dalphac)/eps_c0)**(n*k))) * (alphac+dalphac) - sigc)/dalphac
    else: 
        sigc = 0.1
        dsigc_dalphac = 0.
        
    return sigc,dsigc_dalphac 


def Ds_calc(E_0, ft, fc, eps_c0, x, beta_E, eps_asr, alpha, alphac, sig):  
    sigc = sigc_calc(fc,E_0,eps_c0,alphac)[0]
    Celd = alphac/sigc #+ eps_asr/(beta_E*E_0)
  
    Cs = Celd
    
    
    if sig > 0.:
        sigcr = sigcr_calc(E_0, ft, 0., 0., 0., alpha, 0.)[0]
        Ccr = alpha/sigcr
        Cs = Cs + Ccr

    
    Ds = 1./Cs
    
    return Ds    
    
    