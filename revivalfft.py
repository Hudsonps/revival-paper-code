# -*- coding: utf-8 -*-
from pylab import *
from matplotlib import *
import numpy as np
from math import *
from matplotlib import rc
matplotlib.rcParams.update({'font.size': 17})
from scipy import stats
from scipy import signal

Pi = np.pi
I = 1j

#Amplitude distribution
def coef(n, alpha):
    
    ans = np.sqrt(poisson(n, alpha*alpha))
    return ans
    #return np.sqrt(poisson(n, alpha*alpha))
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
    for n in range(0, 100):
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

def W0_theoretical(alpha, coef, nu):
    W0 = nu
    W0 = 0
    for n in range(0, 100):
        W0 = W0 + Pi*Pi*np.sqrt(nu*nu)*np.sinc(Pi*Pi*nu*nu - n - 1)*coef(n,alpha)*coef(n,alpha)
    return W0

def W1_theoretical(W0theor, nu):
    return 2.*W0theor*np.cos(2.*Pi*Pi*Pi*nu*nu)

def W2_theoretical(W0theor, nu):
    return 2.*W0theor*np.cos(4.*Pi*Pi*Pi*nu*nu)
    

t_0 = -500.
t_f = 500
N = 2**13 #how many slices in the time range
t_ = np.arange(t_0,t_f, (t_f-t_0)/(1.*N)) #vector storing the times to be sampled
deltat = t_[2] - t_[1] #smallest time increment
print deltat
taux_ = reorder(t_) #reorder the vector to make it work with fft properly
#print taux_
alpha = math.sqrt(10) #root of the average photon number
alpha2 = math.sqrt(25)


#time window for the collapse
t1_0 = -10
t1_f = 10
#time window for the different revivals
t2_0 = 5
t2_f = 40
#
t3_0 = 4.4
t3_f = 16
#
t4_0 = 16.
t4_f = 32


ytotalnormal_ = population(alpha, coef, t_)

plt.figure()
plt.xlim(-12, 12)
plt.ylim(-1.3, 2.0)
plot(t_, ytotalnormal_, '-', linewidth = 1.0, color = 'red')
plot(t_, rect(t_, -5, 5), '-', linewidth = 2.0, color = 'coral', label = r'$\Pi_{2T_A}(t)$')
plot(t_, rect(t_, -5, 5) - rect(t_, -2, 2), '--', linewidth = 3.0, color = 'blue', label =  r'$\Pi_{2T_A}(t)$ - $\Pi_{2T_B}(t)$')
plt.xlabel(r'$gt$')
plt.ylabel(r'$W (t) $')
plt.legend(frameon = False)

freq_ = np.fft.fftfreq(N, deltat)
W0theor = W0_theoretical(alpha, coef, freq_)
W1theor = W1_theoretical(W0theor, freq_)
W2theor = W2_theoretical(W0theor, freq_)
T = t_[N-1]-t_[0]
B = N/T
W0time = B*np.fft.ifft(W0theor)
W0time = reorder(W0time)
W1time = B*np.fft.ifft(W1theor)
W1time = reorder(W1time)
W2time = B*np.fft.ifft(W2theor)
W2time = reorder(W2time)


plt.figure()
plt.subplot(2, 1, 1)
plt.xlim(-0, 24)
plt.ylim(-0.6, 1.3)
plot(t_, ytotalnormal_, '-', linewidth = 2.0, color = 'red', label = r'$W(t)$')
#plot(t_, W0time, '--', linewidth = 2.0, color = 'blue', label = r'$W_0(t)$')
#plot(t_, W1time, '--', linewidth = 2.0, color = 'coral', label = r'$W_1(t)$')
plot(t_, W1time + W0time, '--', linewidth = 2.0, color = 'blue', label = r'$W_0 (t) + W_1(t)$')
#plt.xlabel(r'$gt$')
plt.legend(frameon = False,prop={'size':12})


plt.subplot(2, 1, 2)
plt.xlim(-0, 24)
plt.ylim(-0.6, 1.05)
#plot(t_, ytotalnormal_, '-', linewidth = 2.0, color = 'red', label = r'$W(t)$')
plot(t_, W0time, '--', linewidth = 2, color = 'blue', label = r'$W_0(t)$')
plot(t_, W1time, '-.', linewidth = 3.0, color = 'black', label = r'$W_1(t)$')
plot(t_, W2time, ':', linewidth = 3, color = 'magenta', label = r'$W_2(t)$')
#plot(t_, W1time + W0time, '--', linewidth = 2.0, color = 'blue', label = r'$W_0 (t) + W_1(t)$')
plt.xlabel(r'$gt$')
plt.legend(frameon = False,prop={'size':12})



plt.figure()
plot(t_, ytotalnormal_, '-')

#ytotalnormal2_ = population(alpha, coef2, t_)


#we multiply the atomic population by the corresponding rectangle functions to isolate a revival
ytotal_ = population(alpha, coef, taux_)
y1_ = population(alpha, coef, taux_)*rect(taux_, t1_0, t1_f)
y2_ = population(alpha, coef, taux_)*rect(taux_, t2_0, t2_f)
#y3_ = population(alpha, coef, taux_)*rect(taux_, t3_0, t3_f)
#y4_ = population(alpha, coef, taux_)*rect(taux_, t4_0, t4_f)
yft1_ = fft(y1_, N)
yft2_ = fft(y2_, N)
#yft3_ = fft(y3_, t_0, t_f, N)
#yft4_ = fft(y4_, t_0, t_f, N)


freq_ = np.fft.fftfreq(N, deltat) #vector with the frequencies

#plot(freq_, np.abs(yft2_), 'o')
#plot(freq_, np.imag(yft2_), '-')


#inverse transforms
T = t_[N-1]-t_[0]
#yinverse1_ = (N/T)*np.fft.ifft(yft1_)
#yinverse2_ = (N/T)*np.fft.ifft(yft2_)
#yinverse3_ = (N/T)*np.fft.ifft(yft3_)
#yinverse4_ = (N/T)*np.fft.ifft(yft4_)
#plt.xlim(-30,30)
#plot(taux_, yinverse3_ , '')
#plot(taux_,yinverse4_, '-')

yft1m0_ = yft1_
yft1m1_ = yft1_*np.exp(-I*2*Pi*1*((Pi*freq_)**2))
yft1m2_ = yft1_*np.exp(-I*2*Pi*2*((Pi*freq_)**2))
yft1m3_ = yft1_*np.exp(-I*2*Pi*3*((Pi*freq_)**2))
yft1m4_ = yft1_*np.exp(-I*2*Pi*4*((Pi*freq_)**2))

plt.figure()
plot(freq_, fft(ytotal_*rect(taux_, -5, 5), N), color = 'blue')
plot(freq_, fft(ytotal_*rect(taux_, -10, 10), N)/freq_, color = 'red', ls = '--')
plot(freq_, fft(ytotal_*rect(taux_, -20, 20), N), color = 'black')
plot(freq_, fft(ytotal_*rect(taux_, -30, 30), N), color = 'orange')




yinverse1totalm0_ = (N/T)*np.fft.ifft(yft1m0_)
yinverse1totalm1_ = (N/T)*np.fft.ifft(yft1m1_)
yinverse1totalm2_ = (N/T)*np.fft.ifft(yft1m2_)
yinverse1totalm3_ = (N/T)*np.fft.ifft(yft1m3_)
yinverse1totalm4_ = (N/T)*np.fft.ifft(yft1m4_)


#plot(taux_, ytotal_, '-')
#plt.subplot(2, 1, 1)
#plt.title('(a)')
#plt.xlim(0,130)
#plt.ylim(-3, 3)
#plt.xlabel(r'$gt$')
#plt.ylabel(r'$W_m (t)$')
#plot(taux_, yinverse1totalm0_, '-', label=r'$m=0$')
#plot(taux_, yinverse1totalm1_, '-', label=r'$m=1$')
#plot(taux_, yinverse1totalm2_, '-', label=r'$m=2$')
#plot(taux_, yinverse1totalm3_, '-', label=r'$m=3$')
#plot(t_, reorder(yinverse1totalm4_), '-', label=r'$m=4$')
#plt.annotate(r'$W_0 (t)$', xy = (3, 0.9), xytext=(2, 0.65))
#plt.annotate(r'$W_1 (t)$', xy = (25, 0.9), xytext=(25, 0.6))
#plt.annotate(r'$W_2 (t)$', xy = (60, 0.9), xytext=(53, 0.53))
#plt.annotate(r'$W_3 (t)$', xy = (3, 0.9), xytext=(81, 0.44))
#plt.annotate(r'$W_4 (t)$', xy = (3, 0.9), xytext=(107, 0.44))

#plt.axhline
#plt.legend(frameon = False)




    
#plt.subplot(2, 1, 2)
#plt.title('(b)')
#plt.xlim(0,130)
#plt.ylim(-1, 1)
#plt.xlabel(r'$gt$')
#plt.ylabel(r'$W (t)$')
#plot(t_, ytotalnormal_, '-', color = 'red', label = r'$W(t)$')
#plot(t_, reorder(yinverse1total_), '-', label=r'$\sum_{m=0}^4 W_m(t)$')
#plt.tight_layout()



number_ = np.linspace(0,100, 101)
nu_ = np.sqrt(number_+1.)/Pi #those are the special frequencies we need.
zoverf1_ = (1./(2*Pi*Pi))*yft1_/freq_  #fft divided by frequency and 1/pi^2
zoverfinterpolate1_ = interpolate(freq_, zoverf1_, nu_)
zoverf2_ = (1./(2*Pi*Pi))*yft2_/freq_  #fft divided by frequency and 1/pi^2
zoverfinterpolate2_ = interpolate(freq_, zoverf2_, nu_)
#zoverf3_ = (1./(2*Pi*Pi))*yft3_.real/freq_  #fft divided by frequency and 1/pi^2
#zoverfinterpolate3_ = interpolate(freq_, zoverf3_, nu_)
#zoverf4_ = (1./(2*Pi*Pi))*yft4_.real/freq_  #fft divided by frequency and 1/pi^2
#zoverfinterpolate4_ = interpolate(freq_, zoverf4_, nu_)

#using the FFT method to check that the collapses/revival are reproduced
y1prob_ = population2(alpha, zoverfinterpolate1_.real, t_)
y2prob_ = population2(alpha, zoverfinterpolate2_.real, t_)
#y3prob_ = population2(alpha, zoverfinterpolate3_, t_)
#y4prob_ = population2(alpha, zoverfinterpolate4_, t_)

#plot comparing original collapse/revival with one obtained by FFT
#plt.figure()
#plt.xlim(-0,100)
#plt.plot(t_, ytotalnormal_, '-', color = 'red', linewidth = 2.0, label = r'$P_n$')
#plt.plot(t_, y1prob_, '*', color = 'blue', linewidth = 2.0, label = r'$\frac{g^2}{ \pi^2}\frac{\tilde{W_0}(\nu)}{\nu}$')
#xlabel(r'$g t$')
#ylabel(r'$W(t)$')
#plt.legend(frameon=False)

#plt.plot(t_, y3prob_ + y2prob_ + y4prob_, '-', color = 'blue')

print 2.*Pi*np.sqrt(alpha*alpha+1)
plt.figure()
plt.xlim(0,11)
plt.ylim(-1.3,1.3)

#xticks = [2.*Pi*np.sqrt(alpha*alpha+1), 4.*Pi*np.sqrt(alpha*alpha+1), 6.*Pi*np.sqrt(alpha*alpha+1)]
xticks = [Pi, 2*Pi, 3*Pi]
labels = [1, 2, 3] 
plt.xticks(xticks, labels)
plt.plot(t_, ytotalnormal_, '-', color = 'red', linewidth = 2.0, label = r'$| \alpha \rangle$')
#plt.plot(t_, ytotalnormal2_, '-', color = 'blue', linewidth = 2.0, label = r'$|\alpha \rangle$')
#plt.plot(t_, ytotalnormal3_, '-', color = 'green', linewidth = 2.0)
plt.annotate(r"$\langle a^\dagger a \rangle = 20$", xy = (3.14,1.05), xytext = (3.14,1.05), xycoords = "data", textcoords = "data")
plt.xlabel(r"$t/\tau$", fontsize=20)
plt.ylabel(r'$W(t)$', fontsize=20)
plt.legend(frameon=False)



#plt.figure()
#plt.plot(freq_, zoverf1_, '-')
#plt.plot(nu_, zoverfinterpolate1_, 'o')


#plt.figure()
#plt.plot(number_, np.abs(coef(number_, alpha))**2, '--', linewidth = 2.0)
#plt.plot(number_, zoverfinterpolate1_ , 'o', color = 'red', linewidth = 2.0)
#plt.plot(number_, zoverfinterpolate2_ + 2*zoverfinterpolate3_ + zoverfinterpolate4_, 'o', color = 'red')
#plt.plot(number_, poisson(number_, alpha*alpha), 'x')
#plt.plot(number_, zoverfinterpolate2_ , 'h')

plt.figure()
plt.xlim(0,42)
plt.ylim(0,0.2)
x_ = (Pi*freq_)**2 - 1.
plt.plot(number_, np.abs(coef(number_, alpha))**2, 'o', color = 'red', markersize = 4.0, label = r'$P_n$')
plt.plot(x_[:N/2], zoverf1_.real[:N/2] , '--', color = 'blue', linewidth = 2.0, label = r'$\frac{g^2}{ \pi^2}\frac{\tilde{W_0}(\nu)}{\nu}$')
#plt.plot(x_[:N/2], zoverf2_.real[:N/2]/zoverf1_.real[:N/2] , ':', color = 'purple', linewidth = 2.5, label = r'$\frac{g^2}{2 \pi^2}\Re\frac{\tilde{Z_1}(\nu)}{\nu}$')
#plt.annotate(r"$\langle a^\dagger a \rangle = 20$", xy = (31,-0.05), xytext = (31,-0.05), xycoords = "data", textcoords = "data")
plt.xlabel(r'$n$', fontsize=20)
plt.legend(frameon=False)








plt.show()




