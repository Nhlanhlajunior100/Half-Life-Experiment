# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 23:11:38 2021

@author: Nhlanhla Hlengane
"""
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

channels = np.genfromtxt('Co60.tsv', skip_header=22, usecols=(0))
Co_60 = np.genfromtxt('Co60.tsv', skip_header=22, usecols=(1))
Cs_137 = np.genfromtxt('Cs137.tsv', skip_header=22, usecols=(1))
Na_22 = np.genfromtxt('Na22.tsv', skip_header=22, usecols=(1))

alum_spec = np.zeros(1024)
Background = np.zeros(1024)

for i in range(180):
    alum = np.genfromtxt('AluminiumData_'+str(i+1)+'.tsv',skip_header=22,usecols=(1))
    alum_spec+=alum
alum_spec = np.array([float(x)/1800 for x in alum_spec])

for i in range(3):
    Back_10 = np.genfromtxt('Trial1_'+str(i+1)+'.tsv',skip_header=22, usecols=(1))
    Background+=Back_10
Background = np.array([float(x)/30 for x in Background])
alum_spec = alum_spec-Background
Spectra=[Cs_137,Na_22,Na_22,Co_60]
Cent_Guess=[200,200,400,800]
E0=[float(x) for x in Co_60]
Na=[float(x) for x in Na_22]
for i in range(len(Co_60)):
    if channels[i]<400:
        E0[i]=0
    if channels[i]<150:
        Na[i]=0

def Gaus(H,A,sigma,centroid):
    return (A/(sigma*np.sqrt(2*np.pi)))*np.e**(-((H-centroid)**2)/2*sigma**2)
PPopt=[]
PPcov=[]
name=['Caesium-137','Sodium-22','Sodium-22','Cobalt-60']
for i in range(4):
    initialGuess=[1,1,Cent_Guess[i]]
    plt.scatter(channels,Spectra[i])
    #Cobalt60 fix
    if i==3:
        Spectra[i]=E0
    
    #Na_22 Photopeak
    if i==1:
        Spectra[i]=Na
    popt, pcov = curve_fit(Gaus,channels,Spectra[i],initialGuess)
   # plt.plot(channels,Gaus(channels,*popt),'r')
   # plt.title(name[i])
   # plt.xlabel('Channel')
   # plt.ylabel('Count')
   # plt.show()
    PPopt.append(popt)
    PPcov.append(pcov)

for i in range(4):
    print("The peak centroid for %s is at %f"%(name[i],PPopt[i][2]))

#linearized function
def Least(t,m,c):
    return np.log(c)+m*t

    
#Energy calibration    
Known_E=[661.7,1274.5,511,1332]
Cent=[208.728644,403.919507,161.381309,421.651709]
#plt.scatter(Cent,Known_E)
#plt.grid()
#plt.title("Energy Calibration")


m=-173
c=2
p1=[m,c]
name1=['gradient','intercept'] #list of parameter names

popt1,pcov1=curve_fit(Least,Cent,Known_E,p1,absolute_sigma=True)
x =np.linspace(0,500)

y_lin=Least(x,*popt1)
plt.plot(x,y_lin,'r')
plt.show()

#Calibrated Spectrum
M=popt1[0]
C=popt1[1]
print('calibration gradient is '+str(popt1[0])+' and the intercept is '+str(popt1[1]))

#Cal_spec=[popt1[0]*x+popt1[1] for x in alum_spec]
#Sodium22 = [popt[0]*x+popt[1] for x in channels]
Energy = []
for i in range(len(channels)):
    Energy.append(M*channels[i]+C)
plt.step(Energy,alum_spec)
plt.title("Aluminium-28")
plt.xlabel("Energy(keV)")
plt.ylabel("Counts per second")
plt.xlim([0,3300])
#plt.ylim
plt.yscale('log')
plt.grid()
plt.show()

subx=[]
suby=[]
for i in range(len(Energy)):
    if Energy[i]>1400 and Energy[i]<2000:
        subx.append(Energy[i])
        suby.append(alum_spec[i])
plt.step(subx,suby)
plt.title("Aluminium-28")
plt.xlabel("Energy(keV)")
plt.ylabel("Counts per second")
#plt.xlim([0,3300])
#plt.ylim
plt.yscale('log')
plt.grid()
plt.show()
       
################################################3
Background=[]
Back_1 = np.genfromtxt('Trial1_'+str(1)+'.tsv',skip_header=22, usecols=(1))
Back_2 = np.genfromtxt('Trial1_'+str(2)+'.tsv',skip_header=22, usecols=(1))
Back_3 = np.genfromtxt('Trial1_'+str(2)+'.tsv',skip_header=22, usecols=(1))
for i in range(len(Back_1)):
    Background.append((Back_1[i]+Back_2[i]+Back_3[i])/3)
Background=np.array([float(x)/10 for x in Background]) #Background counts per second

CR=[]
region=[]
lenlist = []

for i in range(72):
    region=[]
    x=[]
    al = np.genfromtxt('AluminiumData_'+str(i+1)+'.tsv',skip_header=22,usecols=(1))
    al = np.array([float(x)/10 for x in al])
    al = (al-Background) #subtracting counts per second
    for j in range(len(Energy)):
        if Energy[j]>=1400 and Energy[j]<=1900:
            if al[j]>0:
                region.append(al[j])
                x.append(j)
    CR.append(region)
    lenlist.append(len(region))
    print("Sum of data:",sum(region))

R=1800-1500
Decaylist = [] #region data points
Timelist = []

count=1
for i in range(R):
    a=[]
    for j in range(len(CR)):
        if i<len(CR[j]) and CR[j][i]>0:
            a.append(CR[j][i])  
    if len(a)>0: 
        Decaylist.append(a)


for i in range(len(Decaylist)):
    time = np.linspace(0,1800,len(Decaylist[i]),endpoint=True)
    Timelist.append(time)

Ppopt=[]
Ppcov=[]

m=-173
c=2
p1=[m,c]
U=[]
YLin=[]

for i in range(len(Decaylist)):
    y_lin=np.zeros(len(Decaylist[i]))
    u=np.zeros(len(y_lin))
    for j in range(len(Decaylist[i])):
        if Decaylist[i][j]!=0:
            y_lin[j]=np.log(Decaylist[i][j])
            u[j]=np.sqrt(Decaylist[i][j])/(Decaylist[i][j])
    U.append(u)
    YLin.append(y_lin)
   
   
for i in range(len(Decaylist)):
    if len(YLin[i])<20:
        break
    else:
        popt2,pcov2=curve_fit(Least,Timelist[i],YLin[i],p1,sigma=U[i],absolute_sigma=True)
        Ppopt.append(popt2)
        Ppcov.append(pcov2)
       # plt.scatter(Timelist[i],YLin[i])
        #plt.scatter(Timelist[i],Decaylist[i])
       # plt.plot(Timelist[i],Least(Timelist[i],*popt2),'r')
       # plt.show()
        chisq=sum(((YLin[i]-Least(Timelist[i],*popt2))**2)/U[i])
        dof=len(YLin[i])-len(popt2)
       # print("CHI Squared per degree of freedom:",chisq/dof)
        #print("time",i)
       # print("data points:",len(YLin[i]))
        
        
L=[]
R0=[]
Err=[]
for i in range(len(Ppopt)):
    L.append(Ppopt[i][0])
    R0.append(Ppopt[i][1])
    Err.append(np.sqrt(np.diag(Ppcov[i])))
L=[-1*x for x in L]
T_Half = np.array([np.log(2)/x for x in L])

#intercept=sum(R0)/len(R0)
#print(T_Half)
print("----------------------------------------------------------------")
#print("The average half life of 28AL using this data set is:",datetime.timedelta(seconds=np.mean(T_Half)))
#print("The average half life of 28AL using this data set is:",np.mean(R0))    
#Goodness of fit


Sum=np.zeros(len(CR))
for i in range(len(CR)):
    Sum[i]=sum(CR[i])
Sum=[x for x in Sum]
T=np.linspace(0,720,len(Sum))
SumLin=[np.log(x) for x in Sum]   
unc=[np.sqrt(x)/x for x in Sum]
unc1=[np.sqrt(x) for x in Sum]

popt,pcov=curve_fit(Least,T,SumLin,p1,sigma=unc,absolute_sigma=True)
print("Okay I hope this works! Your Half life is:",np.log(2)/popt[0])
print(np.sqrt(np.diag(pcov)))        
plt.plot(T,popt[1]*np.e**(popt[0]*T), label='28Al Decay model')
plt.errorbar(T, Sum, yerr=unc1,fmt='.', marker='s', ms=2, ecolor='k', mfc='red', mec='red', elinewidth=1, capsize=2, capthick=1, label='28Al Decay')
plt.grid()
plt.title("Aluminium-28 Decay")
plt.legend()
plt.xlabel("Time(seconds)")
plt.ylabel("Counts per second")
plt.xlim([0,730])