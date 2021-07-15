import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import figure
import matplotlib.font_manager
from os import system
from sys import argv
from matplotlib import rc
#import csv
#import pandas as pd
import os
import subprocess
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
#fig = mpl.pyplot.figure(figsize=(7,5))

#fig = mpl.pyplot.figure(figsize=(12,9.5))
#ax1 = fig.add_subplot(111)


n = 10000
fermi = 12.7987

def main(argv):
  
  ymax = 4
  ymin = 0
  #fermi = subprocess.call("grep Fermi qe-scf*.out | awk 'BEGIN{temp=$5}{temp=$5}END{print temp}'", shell=True)
  print("The detected Fermi energy is " + str(fermi) + " eV")
  x = []
  s_up = []
  s_down = []
  p_up = []
  p_down = []
  d_up = []
  d_down = []
  atomssubplot = [[2,3,4,5,6,7,8],[1],[1,2,3,4,5,6,7,8]] 
  #atomssubplot = [[4,9],[1,6],[16],[11]]
  directory = './'
  #bigax = fig.add_subplot(111)
  fig, axs = plt.subplots(len(atomssubplot),sharex=True, sharey=True, gridspec_kw={'hspace': 0.1}) 
  #mpl.pyplot.figure(figsize=(7,5))
  i = 0 
  for layer in atomssubplot: 	#LOOP OVER ALL ATOMS IN GIVEN LAYER
   print("Processing atoms: " + str(layer))
   x = np.zeros(n)
   s_up = np.zeros(n)
   s_down = np.zeros(n)
   p_up = np.zeros(n)
   p_down = np.zeros(n)
   d_up = np.zeros(n)
   d_down = np.zeros(n)
   for atom in layer:		#ANALYSE EACH ATOM 
    prefix = 'scn.pdos.dat.pdos_atm#' + str(atom) + "("
    print(prefix)
    for filename in os.listdir(directory): 
     if filename.startswith(prefix): 			#LOOP OVER EACH ORBITAL FILE
      with open(filename) as fp:		
        print(filename)
        line = fp.readline()
        line = fp.readline()
        j = 0
        x = np.zeros(n)
        yup = np.zeros(n)
        ydown = np.zeros(n)
        while (line and j < n):				#
           line = fp.readline()
           columns = line.split()
           if(len(columns)>1):
            #np.append(x,(float(columns[0]) - fermi))
            x[j] = (float(columns[0]) - fermi)
            yup[j] = (float(columns[1]))
            ydown[j] = (-float(columns[2]))
            j = j + 1
            color = '#000000'
        if filename.endswith('s)'):
              #s_up = (np.add(s_up, yup)).tolist()
              #s_down = (np.add(s_down, ydown)).tolist()
          s_up = np.add(s_up, yup)
          s_down = np.add(s_down, ydown)
        if filename.endswith('p)'):
          p_up = np.add(p_up, yup)
          p_down = np.add(p_down, ydown)
              #p_up = (np.add(p_up, yup)).tolist()
              #p_down = (np.add(p_down, ydown)).tolist()
        if filename.endswith('d)'):
          d_up = np.add(d_up, yup)
          d_down = np.add(d_down, ydown)
            #axs[i].plot(x,yup,c=color)
            #axs[i].plot(x,ydown,c=color)
            #axs[i].plot([0,0],[-10,10],c='#444444', lw = 0.5)
            #axs[i].set_xlim(-10.5,10.5)
            #axs[i].set_ylim(-2.3,2.3)
   #print(str(x[5000]) + " " + str(p_up[5000]))
   total_up = np.add(s_up,p_up)
   total_up = np.add(total_up,d_up)
   total_down = np.add(s_down,p_down)
   total_down = np.add(total_down,d_down)
   axs[i].plot(x,s_up,c='red',label="s")
   axs[i].plot(x,s_down,c='red')
   axs[i].plot(x,p_up,c='limegreen',label="p")
   axs[i].plot(x,p_down,c='limegreen')
   axs[i].plot(x,d_up,c='blue',label="d")
   axs[i].plot(x,d_down,c='blue')
   #axs[i].plot(x,total_up,c='black',label="Total")
   #axs[i].plot(x,total_down,c='black')
   axs[i].plot([0,0],[-10,10],c='#444444', lw = 0.5)
   axs[i].set_xlim(-4,4)
   axs[i].set_ylim(-6,6)
   i = i + 1   
  #titles = ['5','5','5','5','5','5','5','5','5','5','5','5','5','5','5','5','"Bulk" ScN','Interface ScN','Interface Fe','Bulk Fe']
  #titles = ['','','','','','','','','','','','','','','','','','','','','','','','','','']

  titles = ['Sc' + r'$_{0.75}$' + 'N PDOS', 'Fe' + r'$_{0.25}$' + ' PDOS', 'Fe' + r'$_{0.25}$' + 'Sc' + r'$_{0.75}$' + 'N Total DOS']
  for ax in axs:
    #ax.set_yticks(np.arange(-1, 1.1, step=0.5))
    ax.set_yticks(np.arange(-6, 6.1, step=3))
    ax.set_xticks(np.arange(-4, 4.1, step=2))
    ax.label_outer()
    ax.arrow(-6.1, 1, 0, 0.8, head_width=0.2, head_length=0.2, fc='k', ec='k')
    ax.arrow(-6.1, -1, 0, -0.8, head_width=0.2, head_length=0.2, fc='k', ec='k')
  for i in range(len(axs)):
    #axs[i].set_ylabel(titles[i],font='Arial',fontweight='bold',fontsize=12)
    axs[i].set_ylabel(titles[i],font='Arial',fontsize=12)
  #axs[0].set(ylabel=titles[i])
  #axs[int(len(axs)/2)].yaxis.set_label_position("right")
  #yonright = fig.twin()
  axs[len(axs)-1].set_xlabel("E - E" + r'$_f$' + " (eV)",font='Arial',fontsize=15)  
  #axs[len(axs)-1].set_xlabel("E - Ef (eV)",font='Arial',fontweight='bold',fontsize=15)  
  #yonright.set_ylabel("Partial Density of States")
  plt.gcf().subplots_adjust(left=0.12)
  plt.gcf().subplots_adjust(bottom=0.12)
  fig.set_size_inches(8,9)
  #fig.suptitle('Local DOS vs Energy (eV) in Fe/ScN[6]/Fe',y=0.95,font='Arial',fontweight='bold',fontsize=20)
  axs[0].legend(loc='upper right', borderaxespad=0.2)
  plt.rcParams["mathtext.fontset"] = "stix"
  plt.savefig("dos.png")
  plt.show()


#def readfile(filename,

main(argv)   

