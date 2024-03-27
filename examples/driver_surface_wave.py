#This set of codes:
# 1. Determines the root(s) of dispersion equation 
#   for Love waves in a two-layer system by calling the root_newton_raphson function
#   for different frequencies and modes
# 2. Plots figures of solutions, Love wave velocities and wavelengths vs frequencies
#    and of the function in question

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pretty_errors
import os
figpath = "/Users/arvinkarpiah/Desktop/GOPH420/Lab_02/Repository/goph420_w2024_lab02_AK/figures"

from my_python_package.rootfind import root_newton_raphson

# Define frequency
frequency = np.array([0.25, 0.5, 1, 2, 5, 10])

# Define zeta_max
beta1 = 1900
beta2 = 3200
H = 4000
zeta_max = np.sqrt((H ** 2)*((beta1 ** -2)-(beta2 ** -2)))

# Initialize vectors
x_values = np.linspace(0.01, 1.68, num=1000)
root = np.zeros((len(frequency), 3))
num_iter = np.zeros((len(frequency), 3))
zeta_initial = zeta_max

# Define the function g and its derivative
def g(x):
        global rho1
        global rho2
        global beta1
        global beta2
        global H
        rho1 = 1800
        rho2 = 2500
        beta1 = 1900
        beta2 = 3200
        H = 4000
        f = frequency[i]
        return np.tan((2 * np.pi * f) * x) - (rho2 / rho1) * (
                (np.sqrt((H ** 2) * ((beta1 ** -2) - (beta2 ** -2)) - x ** 2)) / x)

def g_prime(x):
    
        f = frequency[i]
        S = H ** 2 * (beta1 ** -2 - beta2 ** -2)
        term1 = 2 * np.pi * f * (1 / np.cos(2 * np.pi * f * x) ** 2)
        term2 = (rho2 / rho1) * (np.sqrt(S-x**2))/x**2
        term3 = (rho2 / rho1) * (1/(np.sqrt(S-x**2)))

        return (term1 + term2 + term3)

# Determine initial guess
for i in range(len(frequency)):
    
        k = 0  # Set mode
        while k < 3:  
            zeta_asymptotes = (0.25 / frequency[i]) * (2 * k + 1)
            if zeta_asymptotes > zeta_max:
                break
            else:
                zeta_initial = zeta_asymptotes - 1e-4
                root[i, k], num_iter[i,k], residual = root_newton_raphson(zeta_initial, g, g_prime)
                k += 1


#Plot function (only used to test plot one freqeuncy at a time)
# plt.figure
# plt.figure(figsize=(10, 6))
# plt.plot(x_values, g(x_values),'o',label='0.5Hz')
# plt.title('Plot of function g')
# plt.xlabel('Time (s)')
# plt.axvline(x=zeta_asymptotes, color='r', linestyle='--',label='asymptote')
# plt.ylim(-10,10)
# plt.ylabel('g')
# plt.legend()
# filename = "figure0.png"
# filepath = os.path.join(figpath, filename)
# plt.savefig(filepath)
                
# Convert the matrix to a pandas DataFrame
df = pd.DataFrame(root)

# Export the DataFrame to a CSV file
df.to_csv('root_table.csv', index=False)               

#Plot solutions against frequencies for different modes
plt.figure
plt.figure(figsize=(10, 6))
plt.plot(frequency, root[:,0],label='mode 0')
plt.plot(frequency[1:], root[1:,1],label='mode 1')
plt.plot(frequency[2:], root[2:,2],label='mode 2')
plt.title('Plot of zeta values for different frequencies')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Time (s)')
plt.legend()
filename = "figure1.png"
filepath = os.path.join(figpath, filename)
#plt.savefig(filepath)


#Plot Love wave velocities against frequencies for different modes
plt.figure 
plt.figure(figsize=(10, 6))
plt.plot(frequency, ((beta1 ** -2) - ((root[:,0]/H) ** 2)) ** -0.5,label='mode 0')
plt.plot(frequency[1:], ((beta1 ** -2) - ((root[1:,1]/H) ** 2)) ** -0.5,label='mode 1')
plt.plot(frequency[2:], ((beta1 ** -2) - ((root[2:,2]/H) ** 2)) ** -0.5,label='mode 2')
plt.title('Plot of Love wave velocities for different frequencies')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Velocity (m/s)')
plt.legend()
filename = "figure2.png"
filepath = os.path.join(figpath, filename)
#plt.savefig(filepath)


#Plot wavelengths against frequencies for different modes
plt.figure
plt.figure(figsize=(10, 6))
plt.plot(frequency, (((beta1 ** -2) - ((root[:,0]/H) ** 2)) ** -0.5)/frequency ,label='mode 0')
plt.plot(frequency[1:], (((beta1 ** -2) - ((root[1:,1]/H) ** 2)) ** -0.5)/frequency[1:],label='mode 1')
plt.plot(frequency[2:], (((beta1 ** -2) - ((root[2:,2]/H) ** 2)) ** -0.5)/frequency[2:],label='mode 2')
plt.title('Plot of wavelengths for different frequencies')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Wavelength (m)')
plt.legend()
filename = "figure3.png"
filepath = os.path.join(figpath, filename)
#plt.savefig(filepath)

#Plot number of iterations against frequencies for different modes
plt.figure
plt.figure(figsize=(10, 6))
plt.plot(frequency, num_iter[:,0],label='mode 0')
plt.plot(frequency[1:], num_iter[1:,1],label='mode 1')
plt.plot(frequency[2:], num_iter[2:,2],label='mode 2')
plt.title('Plot of Number of iterations to converge for different frequencies')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Number of iterations to converge')
plt.legend()
filename = "figure4.png"
filepath = os.path.join(figpath, filename)
#plt.savefig(filepath)
#plt.show()






