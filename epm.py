import numpy as np
from numpy import linalg as LA
import scipy as sci
from scipy import linalg as SLA
import matplotlib as mpl
import matplotlib.pyplot as plt

def main():
    #some constants
    a= 5.431 /0.529; 
    latA=(2*np.pi/a)**2;
    Gmax = 5
    Ecut = 15.5
    print("Analyzing Si with a = " + str(a) + "a.u.")
    print("Ecut = " + str(Ecut) + " Hartrees, Gmax = " + str(Gmax))
    #points to plot between
    X = np.array([0,0,1])
    U = np.array([1/4,1/4,1])
    K = np.array([0,3/4,3/4])
    G = np.array([0,0,0])
    L = (1/2)*np.array([1,1,1])
    kcoeff = 1
    highSymPoints = np.array([L,G,X,U,K,G])
    pathDensity =128
    kGrid = interpolate(highSymPoints,pathDensity)
    gGrid = gPointsGen(Gmax)
    E = []
    #loop over k points and find eigvls of hamiltonian
    for k in kGrid:
        print("Analyzing k = " + str(k/kcoeff))
        kgGrid = kgGridGen(k,gGrid,Ecut)
        H = hamiltonianDiagonals(kgGrid,latA)
        H = hamiltonianOffDiagonals(kgGrid,H)
        eigvls= SLA.eigvalsh(H)
        Ef = 9.5623 + 1.3
        relE = 27.2*eigvls - Ef
        relE.sort()
        E.append(relE)
    printEnergies(E,kGrid,kcoeff)
    #print("Fermi Estimate: " + str(Ef))
    answerQuestions(E)
    plotBands(E,pathDensity,a)

def plotBands(E,pathDensity,a):
    fig = mpl.pyplot.figure(figsize=(7,5))
    ax1 = fig.add_subplot(111)
    size = len(E) - 1
    k = []
    ymax = 10
    for i in range(0,20):
        x = []
        Evals = []
        for j in range(0,size):
            xval = j/pathDensity 
            x.append(xval)
            Evals.append(E[j][i])
        plt.plot(x,Evals, lw = 2)
        #plt.scatter(x,Evals, c=Evals, s = 5, marker='.', vmin=-ymax, vmax=ymax, cmap = 'coolwarm')
    y = [-ymax,ymax]
    kpts = [0,1,2,3,4]
    kptnames = ['L','G','X','U,K','G']
    for kpt in kpts:
        x = []
        x.append(kpt)
        x.append(kpt)
        plt.plot(x,y,c='#111111')
    plt.xticks(kpts,kptnames)
    plt.xticks(kpts,kptnames)
    ax1.set_ylabel("E - Ef (eV)", fontsize=15)
    ax1.set_xlabel("k", fontsize=15)
    ax1.set_ylim(-ymax,ymax)
    ax1.set_xlim(0,4)
    plt.show()

def answerQuestions(Eofk):
    vMax = -20
    vMin = 0
    cMin = 20
    for Eatk in Eofk:
        for occupation in range(0,5):
            if(occupation == 3):
                if(Eatk[occupation] < vMin):
                    vMin = Eatk[occupation]
            if(occupation == 3):
                if(Eatk[occupation] > vMax):
                    vMax = Eatk[occupation]
            if(occupation == 4):
                if(Eatk[occupation] < cMin):
                    cMin = Eatk[occupation]
    print("\nValence Band Width: " + "{:.2f}".format(vMax-vMin) + " eV")
    print("\nBand gap: : " + "{:.2f}".format(cMin - vMax) + " eV")
    #print("\n" + "{:.2f}".format(vMax)+ ", " + "{:.2f}".format(cMin) + ", " + "{:.2f}".format(cMin - vMax) + " eV")



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

def gPointsGen(radius):
    gpoints = []
    for i in range(-radius-1,radius+1):
        for j in range(-radius-1,radius+1):
            for k in range(-radius-1,radius+1):
                if(i**2+j**2+k**2 <= radius**2):
                    if((i%2==1 and j%2==1 and k%2==1) or (i%2==0 and j%2==0 and k%2==0)):
                        g = np.array([i,j,k])
                        gpoints.append(g)
    return gpoints


def kgGridGen(k,ggrid,Ecut):
    kgGrid = []
    for i in ggrid:
        kg = np.add(i,k)
        if((np.dot(kg,kg)) <= Ecut):
            kgGrid.append(kg)
    return kgGrid


def hamiltonianDiagonals(kgGrid,latA):
    diaglist = []
    for kg in kgGrid:
        H = np.dot(kg,kg) * latA / 2 
        diaglist.append(H)
    size = len(diaglist)
    H = np.zeros((size,size))
    for i in range(0,size):
        H[i,i] = diaglist[i]
    return H

def powerMethod(H):
    s = len(H)
    #print(str(H))
    Heff = H
    eigvls = []
    for i in range(0,20):
        v1 = np.random.rand(s)
        v2 = v1*2
        while(norm(v1-v2) > 0.000001):
            v2 = v1
            v1 = np.matmul(Heff,v1) / norm(v2)
        v2 = np.matmul(Heff,v1)
        eigvl = norm(v2)/norm(v1)
        eigvls.append(eigvl)
        print(v1)
        diff = eigvl*np.outer(v1,v1)
        print(eigvl)
        print(diff)
        print(Heff)
        Heff = Heff - diff
    #print(str(Heff))
    return np.array(eigvls)

def norm(vec):
    return LA.norm(vec)

def hamiltonianOffDiagonals(kgGrid,H):
    v3 = -0.21/2
    v8 = 0.04/2
    v11 = 0.08/2
    #v3 = -0.23/2
    #v8 = 0.03/2
    #v11 = 0.08/2
    size = len(kgGrid)
    for i in range(0,size-1):
        for j in range(i+1,size-1):
            g = np.add(kgGrid[i],(-1)*kgGrid[j])
            sg = np.cos((1/4)*np.pi*(g[0]+g[1]+g[2]))
            test = g.dot(g)
            #minimize branches
            v = 0
            test3 = np.abs(test-3)
            test8 = np.abs(test-8)
            test11 = np.abs(test-11)
            #if(test3*test8*test11 < 1):
            if(test3 < 0.01):
                v = v3
            elif(test8 < 0.01):
                v = v8
            elif(test11 < 0.01):
                v = v11
            else:
                v = 10**(-10)
            H[i,j] = sg*v
            H[j,i] = H[i,j]
    return H

def printEnergies(E,kGrid,kcoeff):
    CBM = 20
    VBM = -20
    print("\nEnergies vs k:")
    for i in range(0,len(E)):
        eigvls = E[i]
        print("\nK-point: " + str(kGrid[i]/kcoeff))
        print("Energy (eV):")
        for j in eigvls:
            if(abs(j) < 6):
                if(j > 0 and j < CBM):
                    CBM = j
                if(j < 0 and j > VBM):
                    VBM = j
                print("{:.2f}".format(np.real(j)), end =" ")
        print(" ")
    print("\nValence Max, Conduction Min, band gap: " + "{:.2f}".format(VBM)+ ", " + "{:.2f}".format(CBM) + ", " + "{:.2f}".format(CBM - VBM) + " eV")

def fermiEstimate(lat,N):
    eV=13.6*2
    vol = (lat**3)/4
    return eV*(1/2)*((3*np.pi**2)*N/vol)**(2/3)
main()














def kpointsGen(gridSize):
    #not symmetrizing grid.. easy
    kpoints = []
    for i in np.arange(-0.5,0.5,1/gridSize):
        for j in np.arange(-0.5,0.5,1/gridSize):
            for k in np.arange(-0.5,0.5,1/gridSize):
                kpoint = np.array([i,j,k])
                kpoints.append(kpoint)
    return kpoints




