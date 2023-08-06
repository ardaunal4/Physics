import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad 
import seaborn as sns
sns.set()

c = 299792458 # m/s
KE = 2 # MeV
v = c*np.sqrt(1-1/(KE/938.2720813+1)**2) # m/s
KP = 15 # MV/m according to kilpatrick
beta = v/c # need to be evaluated
Lcell = beta*1.5 # lambda =1.5

Lgap = np.linspace(start = 0.0001, stop = 0.2, num = 2000) # m

max_TT = 0
max_Eavg = 0
max_Vo = 0
max_gap_length = 0
max_sigma = 0
optimized_Eo = 0
optimized_z = 0

for gap_length in Lgap:

    z = np.linspace(start = -gap_length/2, stop = gap_length/2, num = 2000)

    for sigma in z:

        Eo = np.exp(-z**2/(2*abs(sigma)**2))/(abs(sigma)*np.sqrt(2*np.pi)) # Electric field dist.

        if (max(Eo)) <= KP:

            Vo,_ = quad(lambda z: np.exp(-z**2/(2*abs(sigma)**2))/(abs(sigma)*np.sqrt(2*np.pi)), -gap_length/2, gap_length/2)
            Eavg = Vo/Lcell
            res,_ = quad(lambda z: np.cos(2*np.pi*z/(Lcell*360))*np.exp(-z**2/(2*abs(sigma)**2))/(abs(sigma)*np.sqrt(2*np.pi)), -gap_length/2, gap_length/2)
            TT = res/Vo

            if max_TT < TT:

                optimized_Eo = Eo
                max_TT = TT
                max_Eavg = Eavg
                max_Vo = Vo
                max_gap_length = gap_length
                max_sigma = abs(sigma)
                optimized_z = z

power_per_cell = 1*max_Eavg*max_TT*Lcell*np.cos(2*np.pi*85/360)

print("max_Eo: ",  max(optimized_Eo))
print("Vo: ",  max_Vo)
print("Eavg: ",  max_Eavg)
print("TT: ",  max_TT)
print("Gap Length: ", gap_length)
print("Sigma: ", max_sigma)

plt.figure(figsize=(10, 8))
plt.plot(optimized_z, optimized_Eo)
plt.xlabel("Gap Length(meters)")
plt.ylabel("Electric Field(*10^6)")
plt.grid("on")
plt.show()
