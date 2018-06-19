# -*- coding: utf-8 -*-
from pylab import *
from matplotlib import *
import numpy as np
from math import *
from matplotlib import rc
matplotlib.rcParams.update({'font.size': 17})
import scipy
import numpy.matlib

Pi = np.pi
I = 1j

def Smatrix(t):
    l = len(t)
    dt = t[1] - t[0]
    tmatrix = numpy.matlib.repmat(t, l, 1)
    S = dt*np.cos(Pi/4. - ((tmatrix - tmatrix.T)**2)/(2.*Pi))*np.sqrt(2)/Pi + np.identity(l)
    return S

#Amplitude distribution
def coef(n, alpha):
    return np.sqrt(poisson(n, alpha*alpha))
    #return np.exp(-0.02*np.abs(n-15)) + np.exp(-0.02*np.abs(n-15))*(np.exp(180*I*Pi/180)**(n))
    #beta = alpha
    #return (np.sqrt(poisson(n, alpha*alpha)) + np.sqrt(poisson(n, beta*beta))*(np.exp(180*I*Pi/180)**(n)))/np.sqrt(2)
         #+ np.sqrt(poisson(n, beta*beta))*(np.exp(100*I*Pi/180)**(n))+ np.sqrt(poisson(n, beta*beta))*(np.exp(180*I*Pi/180)**(n)))
    #return np.exp(-alpha*np.abs(n))
    #return np.sqrt((1./(alpha+1.))*(((alpha)/(alpha+1.))**n))
    #return 1/np.sqrt(n+1)
    #return (1/np.sqrt(np.cosh(alpha)))*((-np.tanh(alpha))**n)*np.sqrt(1.*(math.factorial(2*n)))/(1.*math.factorial(n)*(2**n))

#Second amplitude distribution    
def coef2(n, alpha):
    return np.sqrt(poisson(n, alpha*alpha))
    
#this function takes the function coef as its argument and returns the corresponding atomic population 
def population(alpha, coef, t_):
    zn_ = t_
    zn_ = 0
    for n in range(0, 1000):
        zn_ = zn_ + np.exp(I*2.*math.sqrt((n+1))*t_)*(np.abs(coef(1.*n,alpha))**2)
        #zn_ = zn_ + np.cos(2.*math.sqrt((n+1))*t_)*(np.abs(coef(1.*n,alpha))**2)
        #zn_ = zn_ + np.cos(2.*n*t_)*(np.abs(coef(1.*n,alpha))**2)
    #return zn_*np.exp(-I*1*t_)
    return zn_

#this function takes a probability vector vector instead of a function
def population2(alpha, prob_, t_):
    zn_ = t_
    zn_ = 0
    N = prob_.size
    for n in range(0, N):
        zn_ = zn_ + np.cos(2.*math.sqrt((n+1))*t_)*prob_[n]
    return zn_

#rectangular function, 1 between x_i and x_f, 0 otherwise.
def rect(x_, x_i, x_f):
    n = x_.size
    ans_ = np.zeros(n)
    #print ans_.size
    for i in range(0, n):
        if(x_[i] > x_i and x_[i] < x_f):
             ans_[i] = 1
    return ans_

#find the index of a vector corresponding to the nearest element of value
def find_nearest_index(array_,value):
    idx = (np.abs(array_-value)).argmin()
    return idx

#given a vector y_, interpolates to get y[x]
def interpolate(t_, y_, x_):
    N = x_.size
    ans_ = np.zeros(N)
    for i in range(0,N,1):
        idx = find_nearest_index(t_, x_[i])
    #we need to check if this point is to the right or to the left
        if((t_[idx] - x_[i]) > 0):
            coef = (y_[idx] - y_[idx-1])/(t_[idx]-t_[idx-1])
            ans_[i] = y_[idx-1] + coef*(x_[i] - t_[idx-1])
        if((t_[idx] - x_[i]) < 0):
            coef = (y_[idx+1] - y_[idx])/(t_[idx+1]-t_[idx])
            ans_[i] = y_[idx] + coef*(x_[i] - t_[idx])
    return ans_
          


def fft(y_, N):
    #total time interval
    T = t_[N-1]-t_[0]
    #smallest time step
    #Time to define the normalization factor.
    A = T/N
    #Taking the Fourier transform and already including the normalization.
    yft_ = A*np.fft.fft(y_)
    #Defining the vector with the corresponding frequencies
    return yft_

#the reordering is just to avoid a phase factor
def reorder(t_):
    N = t_.size
    middle = N/2
    if(t_[middle] != 0):
        print("You are using the function reorder incorrectly!")
    auxplus = t_[middle:N]
    auxminus = t_[0:middle]
    aux = np.concatenate((auxplus, auxminus))
    if(aux.size != N):
        print("You are using the function reorder incorrectly!")
    return aux

def poisson(n, mean):
    return stats.distributions.poisson.pmf(n ,mean)

t_0total = -30
t_ftotal = 30
t_0 = -10
t_f = 10
N = 2**11 #how many slices in the time range
t = np.arange(t_0total,t_ftotal, (t_ftotal-t_0total)/(1.*N)) #vector storing the times to be sampled
deltat = t[2] - t[1] #smallest time increment
print deltat
taux = reorder(t) #reorder the vector to make it work with fft properly
#print taux_
alpha = math.sqrt(1) #root of the average photon number

tau1 = 2*Pi*np.sqrt(alpha*alpha + 1.)
print tau1

W = population(alpha, coef, t)
S = Smatrix(t)

#print "The rank is ", scipy.linalg.matrix_rank(S)

Winitial = W*rect(t, t_0, t_f)
W1 = np.dot(S, Winitial)

plt.figure()
#plot(t, W, "-", linewidth = 2.0)
plot(t, Winitial, "--", color = "red")
#plot(t, W1, "--", color = "orange")
print scipy.linalg.det(S)
#W0 = numpy.linalg.solve(S, Winitial)
#plot(t, W0, "-", color = "purple")

eigenval, eigenveg = scipy.linalg.eig(S)

eig_vals_sorted = np.sort(eigenval)
eig_vecs_sorted = eigenveg[:, eigenval.argsort()]

plt.figure()
plot(eig_vals_sorted, '-')
#plot(eigenval, '-')

plt.figure()
plot(t, eig_vecs_sorted[90], '-')
print eig_vals_sorted[90]
plt.show()
