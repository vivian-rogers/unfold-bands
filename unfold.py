import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import figure

#Creates an array that contains all possible dot products of one atom's kpoints with another (or the same) atom's position
def kToPW(Kvals, Rvals, nx, ny):
    nkpts = len(Kvals)
    PWs = np.zeros((nx * ny, nkpts))
    for iK in range(0, nkpts):
        for iR in range(0, nx * ny):
            k = Kvals[iK]
            r = Rvals[iR]
            PWs[iR,iK] = np.exp((0 + (-1)j) * np.sum(k*r))
    return PWs

def interpolate(highSymPoints,density):
    kpts = []
    U = np.array([1/4,1/4,1])
    K = np.array([0,3/4,3/4])
    for i in range(0,len(highSymPoints)-1):
        begin = highSymPoints[i]
        end = highSymPoints[i+1]
        if(~((np.dot(begin-U,begin-U) < 0.001) and (np.dot(end-K,end-K) < 0.001))):
          for j in range(0,density):
            k = (j/density)*end + (1 - (j/density))*begin
            kpts.append(k)
        #else:
            #kpts.append(K*(1-(1/density)))
    kpts.append(highSymPoints[len(highSymPoints)-1])
    return kpts

def bandpath(nx,ny,Estates,Evals,npoints,nE,Emin, Emax, lamb, theta, SLBZ):
    a = 1
    
    #List of xy pairs to target on path
    if(SLBZ == false):
        Ktargets = ((2/3)*pi)*np.array([[0, 0], [0, 1], [np.tan(pi/6), 1], [0, 0]])
    else:
        n1 = lamb*10**9
        n2 = lamb*10**9
        a1 = np.array([1, 0])
        a2 = np.array([1/2, sqrt(3)/2])
        b1 = n1*a1 + n2*a2
        b2 = -n2*a1 + (n1+n2)*a2
        Ktargets = []
    
    #Generates path from destination points
    Kvals = interpolate(Ktargets,npoints)
    nKvals = len(Kvals)
    #Rvals = a*generateRealSpaceBasis(nx,ny) #gives (x,y) pos per index. Need this function????
    bz1PWs = kToPW(Kvals,Rvals,nx,ny) #PWs for 1st brillouin zone
    bz2PWs = kToPW(Kvals + [0, (4/3)*pi],Rvals,nx,ny) #PWs for 2nd brillouin zone
    DOS = zeros(nKvals, nE)
    
    #move along the xpath
    for i in range len(nKvals):
        #env = envelope(Rvals,Xvals(i,:)) %envelope to project wfc onto site
        #start dissecting the hamiltonian
        PW1 = bz1PWs[:,i]
        PW2 = bz2PWs[:,i]
        for j in range in len(Estates):
            psi = Estates[:,j]
            #performs |<psi|planewave(k)>|^2 == <PW|psi><psi|PW> for changing basis, adds 1st and 2nd BZ
            #overlap = abs((psi')*PW2)^2 + abs((psi')*PW1)^2 ??
            
            #get DOS index and bin the delta(E)*overlap
            iDOS = round(nE*(Evals[j,j] - Emin)/(Emax - Emin))
            if (iDOS > 0 and iDOS <= nE):
                DOS[i,iDOS] += overlap
                
    labels = ['\Gamma','M','K','\Gamma']