# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 01:07:53 2020

@author: simkon
"""
import numpy as np

def sectionLayers(crossSection):
    sectionTypeNumber = int(crossSection[1])
    if sectionTypeNumber == 1:
        h = crossSection[4]
        b = crossSection[5]
        Nsec = int(crossSection[2])
        
        zc_array = np.linspace(-h/2,h/2,2*Nsec+1)
        Ac_array = b*h/Nsec*np.ones(Nsec)
        
        Ns_layers = int(crossSection[3])
        As_array = crossSection[8:8+Ns_layers]
        zs_array = crossSection[8+Ns_layers: 8+2*Ns_layers]
    
    elif sectionTypeNumber == 2:
        hw = crossSection[4]
        hf = crossSection[5]
        bw = crossSection[6]
        bf = crossSection[7]
        Nsec = int(crossSection[2])
        h = hw + hf
        
        dz = int(h/Nsec)
        
        Nsec_f = max(int(hf/dz),1)
        Nsec_w = Nsec - Nsec_f
        dz_f = hf/Nsec_f
        dz_w = hw/Nsec_w
        dA_f = dz_f*bf
        dA_w = dz_w*bw
        
        zc_array = np.zeros(2*Nsec+1)
        for i in range(2*Nsec_f+1):
            zc_array[i] = -h/2 + i*dz_f/2.
        for i in range(1,2*Nsec_w+1):
            zc_array[2*Nsec_f + i] = -h/2 + hf +  i*dz_w/2.
        
        Ac_array = np.ones(Nsec)
        Ac_array[0:Nsec_f] = dA_f*Ac_array[0:Nsec_f]
        Ac_array[Nsec_f:] = dA_w*Ac_array[Nsec_f:]
        
        Ns_layers = int(crossSection[3])
        As_array = crossSection[8:8+Ns_layers]
        zs_array = crossSection[8+Ns_layers: 8+2*Ns_layers]
    elif sectionTypeNumber == 3:
        R = crossSection[4]
        Nsec = int(crossSection[2])
        
        zc_array = np.linspace(-R,R,Nsec+1)
        h = zc_array[1]-zc_array[0]
        b_c = 2.*np.sqrt(R**2. - zc_array**2.)
        Ac_array = np.zeros(Nsec)
        for i in range(Nsec):    
            Ac_array[i] = 0.5*(b_c[i]+b_c[i+1])*h
        
        Ns_layers = int(crossSection[3])
        As_array = crossSection[8:8+Ns_layers]
        zs_array = crossSection[8+Ns_layers: 8+2*Ns_layers]
        
        
    return {'zc':zc_array, 'Ac': Ac_array, 'zs':zs_array, 'As': As_array}