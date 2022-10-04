#-COLLECTION OF FUNCTIONS USED IN THE PROGRAM 
import numpy as np


def Nmatrix(x,z,L):
    N = np.zeros((2,7), dtype=np.float)
    Nwx = np.array([-6*x/L**2 + 6*x**2/L**3, 1-4*x/L + 3*x**2/L**2, 6*x/L**2 - 6*x**2/L**3, -2*x/L + 3*x**2/L**3])
    Nw = np.array([1-3*x**2/L**2+2*x**3/L**3, x-2*x**2/L+x**3/L**2, 3*x**2/L**2-2*x**3/L**3, -x**2/L+x**3/L**2])
    Nu = [1-x/L, x/L]
    Nui = -4/L**2 * x*(x - L)
    N[0,0:2] = Nu
    N[0,2:6] = -z*Nwx
    N[0,6] = Nui
    N[1,2:6] = Nw
    return N

def Bmatrix(x,z,L):
    
    Nux = np.array([-1/L, 1/L])
    Nwxx_x= np.array([-6./L**2+12.*x/L**3., -4./L + 6.*x/L**2, 6./L**2-12*x/L**3, -2./L+6*x/L**2])
    Nuix_x = -8*x/L**2 + 4/L

    B = np.append(Nux, -z*Nwxx_x)
    B = np.append(B,Nuix_x)
    
    return B

def strainCalc_concrete(crossSection,elCoord,De):
    L = np.sqrt( (elCoord[2]-elCoord[0])**2. + (elCoord[3]-elCoord[1])**2.)
    x_array= (L/2.-3.**0.5/6*L, L/2.+3.**0.5/6*L)
    c = (elCoord[2]-elCoord[0])/L
    s = -(elCoord[3]-elCoord[1])/L
    T = np.array([[c , 0., -s , 0., 0., 0., 0], \
                  [0., c , 0., 0., -s, 0.,  0.], \
                  [s , 0., c , 0., 0., 0.,  0.], \
                  [0., 0., 0., 1., 0., 0.,  0.], \
                  [0., s , 0., 0., c , 0.,  0.], \
                  [0., 0., 0., 0., 0., 1.,  0.], \
                  [0., 0., 0., 0., 0., 0.,  1.]])    
    zc_array = crossSection['zc']
    
    de = np.matmul(T,De)

    eps_allIntPt = np.zeros((2,len(zc_array)), dtype=float)
    for i in range(2):
        x = x_array[i]
        for j in range(len(zc_array)):
            zc = zc_array[j]
            B = Bmatrix(x,zc,L)
            eps_allIntPt[i][j] = np.dot(B, de)
    return eps_allIntPt

def strainCalc_steel(crossSection,elCoord,De):
    L = np.sqrt( (elCoord[2]-elCoord[0])**2. + (elCoord[3]-elCoord[1])**2.)
    c = (elCoord[2]-elCoord[0])/L
    s = -(elCoord[3]-elCoord[1])/L
    x_array= (L/2.-3.**0.5/6*L, L/2.+3.**0.5/6*L)
    T = np.array([[c , 0., -s , 0., 0., 0., 0], \
                  [0., c , 0., 0., -s, 0.,  0.], \
                  [s , 0., c , 0., 0., 0.,  0.], \
                  [0., 0., 0., 1., 0., 0.,  0.], \
                  [0., s , 0., 0., c , 0.,  0.], \
                  [0., 0., 0., 0., 0., 1.,  0.], \
                  [0., 0., 0., 0., 0., 0.,  1.]])
    zs_array = crossSection['zs']

    de = np.matmul(T,De)
    eps_allIntPt = np.zeros((2,len(zs_array)), dtype=float)
    for i in range(2):
        x = x_array[i]
        for j in range(len(zs_array)):
            zs = zs_array[j]
            B = Bmatrix(x,zs,L)
            eps_allIntPt[i][j] = np.dot(B, de)
    return eps_allIntPt

#INTERNAL FORCES 
def elInternalForce(crossSection,elCoord, sig_c, sig_s): #shape(sig_c)= (2, number of layers over the height)
    L = np.sqrt( (elCoord[2]-elCoord[0])**2. + (elCoord[3]-elCoord[1])**2.)
    c = (elCoord[2]-elCoord[0])/L
    s = -(elCoord[3]-elCoord[1])/L
    T = np.array([[c , 0., -s , 0., 0., 0., 0], \
              [0., c , 0., 0., -s, 0.,  0.], \
              [s , 0., c , 0., 0., 0.,  0.], \
              [0., 0., 0., 1., 0., 0.,  0.], \
              [0., s , 0., 0., c , 0.,  0.], \
              [0., 0., 0., 0., 0., 1.,  0.], \
              [0., 0., 0., 0., 0., 0.,  1.]])
    x_array= (L/2.-3.**0.5/6*L, L/2.+3.**0.5/6*L)
    zc_array = crossSection['zc']
    Ac_array = crossSection['Ac']
    zs_array = crossSection['zs']
    As_array = crossSection['As']
    
        
    rel_int = np.zeros(7, dtype=np.float64)
    
    #Integration over length with Gauss
    for i in range(2):
        x = x_array[i]
        #Integration over cross section
        for j in range(len(Ac_array)):
            MatMul_1 = Ac_array[j]/6.* (Bmatrix(x,zc_array[2*j],L)*sig_c[i,2*j] + 4.*Bmatrix(x,zc_array[2*j+1],L)*sig_c[i,2*j+1] + Bmatrix(x,zc_array[2*j+2],L)*sig_c[i,2*j+2] )
            rel_int = np.add(rel_int, MatMul_1)
        for k in range(len(zs_array)):
            zs = zs_array[k]
            As = As_array[k]
            B = Bmatrix(x,zs,L)
            MatMul_1 = B * sig_s[i,k] * As
            rel_int = np.add(rel_int, MatMul_1)
    rel_int = rel_int*L/2.
    Rel_int = np.matmul(np.transpose(T), rel_int)
    return Rel_int     

def elStiffMat(crossSection,elCoord,stiff_c,stiff_s):
    L = np.sqrt( (elCoord[2]-elCoord[0])**2. + (elCoord[3]-elCoord[1])**2.)
    c = (elCoord[2]-elCoord[0])/L
    s = -(elCoord[3]-elCoord[1])/L
    T = np.array([[c , 0., -s , 0., 0., 0., 0], \
                  [0., c , 0., 0., -s, 0.,  0.], \
                  [s , 0., c , 0., 0., 0.,  0.], \
                  [0., 0., 0., 1., 0., 0.,  0.], \
                  [0., s , 0., 0., c , 0.,  0.], \
                  [0., 0., 0., 0., 0., 1.,  0.], \
                  [0., 0., 0., 0., 0., 0.,  1.]])
                  
    x_array= (L/2.-3.**0.5/6*L, L/2.+3.**0.5/6*L)
    zc_array = crossSection['zc']
    Ac_array = crossSection['Ac']
    zs_array = crossSection['zs']
    As_array = crossSection['As']
    
    ke = np.zeros((7,7), dtype=float)
    #Integration over length with Gauss
    for i in range(2):
        x = x_array[i]
        #Integration over cross section
        for j in range(len(Ac_array)):
            Ac = Ac_array[j]
            B1 = np.array([Bmatrix(x,zc_array[2*j],L)])
            B2 = np.array([Bmatrix(x,zc_array[2*j+1],L)])
            B3 = np.array([Bmatrix(x,zc_array[2*j+2],L)])
            D1 = stiff_c[i,2*j]
            D2 = stiff_c[i,2*j+1]
            D3 = stiff_c[i,2*j+2]
            MatMul_1 = Ac/6.* (np.matmul(np.transpose(B1),B1)*D1 + 4.*np.matmul(np.transpose(B2),B2)*D2 + np.matmul(np.transpose(B3),B3)*D3 )
            ke = np.add(ke, MatMul_1)
        for k in range(len(zs_array)):
            zs = zs_array[k]
            As = As_array[k]
            B = np.array([Bmatrix(x,zs,L)])
            MatMul_1 = stiff_s[i,k] * B * As
            MatMul_2 = np.matmul(np.transpose(B), MatMul_1)
            ke = np.add(ke, MatMul_2)
    
    ke = ke*L/2.
    Ke = np.matmul(np.matmul(np.transpose(T), ke), T)
    return Ke    


def elExternalForce(elCoord, load):
    L = np.sqrt( (elCoord[2]-elCoord[0])**2. + (elCoord[3]-elCoord[1])**2.)
    x_array= (L/2.-3.**0.5/6*L, L/2.+3.**0.5/6*L)
    c = (elCoord[2]-elCoord[0])/L
    s = -(elCoord[3]-elCoord[1])/L
    T = np.array([[c , 0., -s , 0., 0., 0., 0], \
                  [0., c , 0., 0., -s, 0.,  0.], \
                  [s , 0., c , 0., 0., 0.,  0.], \
                  [0., 0., 0., 1., 0., 0.,  0.], \
                  [0., s , 0., 0., c , 0.,  0.], \
                  [0., 0., 0., 0., 0., 1.,  0.], \
                  [0., 0., 0., 0., 0., 0.,  1.]])    
    Tq  = np.array([[c,-s],
                    [s, c]])
    Qe = load
    qe = np.matmul(Tq,Qe)
    
    #Integration over length with Gauss quadrature
    rel_ext = np.zeros(7, dtype=float)
    
    for i in range(2):
        x = x_array[i]
        z = -715.
        Nmat = Nmatrix(x,z,L)
        MatMul_1 = np.matmul(np.transpose(Nmat), qe)
        rel_ext = np.add(rel_ext, MatMul_1)
    rel_ext = rel_ext*L/2.
    
    Rel_ext = np.matmul(np.transpose(T), rel_ext)
        
    return Rel_ext


def w(x, L, dw):
    return np.dot(np.array([1-3*x**2/L**2+2*x**3/L**3, x-2*x**2/L+x**3/L**2, 3*x**2/L**2-2*x**3/L**3, -x**2/L+x**3/L**2]), dw)

def u(x,L, du):
    return np.dot(np.array([1-x/L, x/L]), du)

def U(u, w, theta):
    return w * np.sin(theta) + u * np.cos(theta)
def V(u,w,theta):
    return -w * np.cos(theta) + u * np.sin(theta)
    
def T(theta):
    s = np.sin(theta)
    c = np.cos(theta)

    return np.array([[c , 0., -s , 0., 0., 0.], \
                     [0., c , 0., 0., -s, 0.], \
                     [s , 0., c , 0., 0., 0.], \
                     [0., 0., 0., 1., 0., 0.], \
                     [0., s , 0., 0., c , 0.], \
                     [0., 0., 0., 0., 0., 1.]])
    
    
def sectionForces(crossSection, sig_c, sig_s):
    zc_array = crossSection['zc']
    Ac_array = crossSection['Ac']
    zs_array = crossSection['zs']
    As_array = crossSection['As']
    
        
    M = 0.
    N = 0.

    for j in range(len(zc_array)-1):
        Ac = Ac_array[j]
        N = N + (sig_c[j]+sig_c[j+1])/2. * Ac

        
    for k in range(len(zs_array)):
        As = As_array[k]
        N = N + sig_s[k] * As

    sig_avg = N/(np.sum(Ac_array)+np.sum(As_array))
    
    sig_c_mod = sig_c - sig_avg
    sig_s_mod = sig_s - sig_avg
    
    for j in range(len(zc_array)-1):
        Ac = Ac_array[j]
        M = M + (sig_c_mod[j]*zc_array[j] + sig_c_mod[j+1]*zc_array[j+1])/2. * Ac

    
    for k in range(len(zs_array)):
        zs = zs_array[k]
        As = As_array[k]
        M = M + sig_s_mod[k] * As*zs

    
    return N,M


def epM(L, Mx1, Mx2,x):
    x1= (L/2.-3.**0.5/6*L)
    x2= (L/2.+3.**0.5/6*L)
    a = (Mx2-Mx1)/(x2-x1)
    b = Mx1-a*x1
    
    return a*x+b
    
def epN(L, Nx1, Nx2,x):
    x1= (L/2.-3.**0.5/6*L)
    x2= (L/2.+3.**0.5/6*L)
    a = (Nx2-Nx1)/(x2-x1)
    b = Nx1-a*x1
    
    return a*x+b 