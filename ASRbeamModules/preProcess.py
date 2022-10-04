# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 17:04:30 2020

@author: simkon
"""

#-PREPROCESSOR
import numpy as np
#import functions as fu
import ASRbeamModules.crossSectionCalc as csc

def readin(folder=""):
    """Read the input files"""
    nodes = np.loadtxt(folder + 'nodes.txt', ndmin=2)
    mat_concrete = np.loadtxt(folder + 'mater_concrete.txt', ndmin=2)
    mat_steel = np.loadtxt(folder + 'mater_steel.txt', ndmin=2)
    elements = np.loadtxt(folder + 'eles.txt', ndmin=2, dtype=np.int)
    crossSections = np.loadtxt(folder + 'crossSections.txt', ndmin=2)
    loads1 = np.loadtxt(folder + 'loads_step1.txt', ndmin=2)
    loads2 = np.loadtxt(folder + 'loads_step2.txt', ndmin=2)
    analysis = np.loadtxt(folder + 'analysis.txt', ndmin=2)
    predef1 = np.loadtxt(folder + 'predef_step1.txt', ndmin=2)
    predef2 = np.loadtxt(folder + 'predef_step2.txt', ndmin=2)

    return nodes, elements, crossSections, mat_concrete, mat_steel, loads1, loads2, analysis, predef1, predef2

def connectivityMat(nodes, elements):
    Nnodes = np.shape(nodes)[0]
    Nel = np.shape(elements)[0]
    
    #Ndofs = 3*Nnodes + Nel
    
    C = np.zeros((Nel, 7), dtype=int)
    for e in range(Nel):
        nodes_ext_el = elements[e][3:]
        
        C[e][0] = nodes_ext_el[0]*3
        C[e][1] = nodes_ext_el[1]*3
        
        C[e][2] = nodes_ext_el[0]*3+1
        C[e][3] = nodes_ext_el[0]*3+2
        C[e][4] = nodes_ext_el[1]*3+1   
        C[e][5] = nodes_ext_el[1]*3+2 
        C[e][6] = 3*Nnodes  + e

    return C


def boundaryCondMat(nodes):
    #--BOUNDARY CONDITIONS
    BC = np.array([], dtype=int)
    
    for i in range(np.shape(nodes)[0]):
        if nodes[i][3] == -1:
            dof = int(nodes[i][0])*3
            BC = np.append(BC,dof)
        if nodes[i][4] == -1:
            dof = int(nodes[i][0])*3+1
            BC = np.append(BC,dof)
        if nodes[i][5] == -1:
            dof = int(nodes[i][0])*3+2
            BC = np.append(BC,dof)
    BC = np.sort(BC)
    return BC

def crossSectionMat(crossSections):
    N = np.shape(crossSections)[0]
    cs_array = np.array([])
    for i in range(N):
        cs_array = np.append(cs_array, csc.sectionLayers(crossSections[i]))
    
    return cs_array

def elCoordMat(nodes,elements):
    Nel = np.shape(elements)[0]
    elCoord = np.zeros((Nel,4), dtype = float)
    for e in range(Nel):
        elCoord[e][0] = nodes[elements[e][3]][1]
        elCoord[e][1] = nodes[elements[e][3]][2]
        elCoord[e][2] = nodes[elements[e][4]][1]
        elCoord[e][3] = nodes[elements[e][4]][2]
        
    return elCoord

 