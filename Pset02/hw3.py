import numpy as np
import matplotlib.pyplot as plt

#Problem 1B
q = 1.60217663 * 10**(-19)
no1 = 10**(13)
no2 = 10**(18)
no3 = 10**(13)
r1 = 0.25*10**(-3)
r2 = 0.5*10**(-3)
r3 = 0.25*10**(-3)
mu1 = 1410
mu2 = 390
mu3 = 1410
r_sh = 1/(q*(no1*r1*mu1+no2*r2*mu2 + no3*r3*mu3))
print("r_sh", r_sh)

#Problem 1C
r_ab = 71*r_sh
print("r_ab", r_ab)

#problem 2B
nd = 10**(1.12/(2*60*10**(-3)))*10**10
print("Nd", '{:.2e}'.format(nd))
rho = 1/(q*nd*125)
print("resistivity", '{:.2e}'.format(rho))

print("Nd", '{:.2e}'.format(nd))
rho = 1/(q*nd*10)
print("resistivity", '{:.2e}'.format(rho))


#Problem 2e
nd = 10**19
na = 10**13
ni = 10**10
epsilon_s = 11.7*8.85*10**(-14)
phi_b = 60*np.log10(na*nd/(ni**2))
print("phib", phi_b)
x_do = np.sqrt(2*epsilon_s/q*phi_b*((na+nd)/(na*nd)))
print("width min", 2*x_do)

#Problem 4b
na = 10**15
nd = 10**21
phi_b = 10**(-3)*(60*np.log10(na*nd/(ni**2)))
print("phib", phi_b)
frac = (4*(10**(-8))*q*na)/(2*epsilon_s)
print("frac", frac)

vt = frac-phi_b
print("vt", vt)

#Problem 6A
L = 2.5e-9
na = 1e16
nd = 1e19
f = 5.2e9
epsilon_0 = 1.04e-12
C = 1/(L*(f*2*np.pi)**2)
phi_b = 60e-3*np.log10(na*nd/(ni**2))
d = np.sqrt(2*epsilon_s*phi_b/q*(na+nd)/(na*nd))
A = C*d/epsilon_0
print("capacitance:", '{:.2e}'.format(C))
print("area", '{:.2e}'.format(A))


#Problem 6b
frequency_list = [5260e6,5280e6,5300e6, 5320e6, 
                  5500e6, 5520e6, 5540e6, 5560e6,5580e6, 
                  5680e6, 5700e6]
voltages = np.linspace(-0.45,0,1000)
frequencies = []
for v in voltages:
    C_prime = C/np.sqrt(1-v/phi_b)
    f = 1/(2*np.pi*np.sqrt(L*C_prime))
    frequencies.append(f)


plt.plot(voltages, frequencies)
plt.xlabel('Voltage (V)')
plt.ylabel('Frequency (Hz)')
plt.title('Voltage vs Frequency')
min_freq = min(frequencies)
max_freq = max(frequencies)
step = round((max_freq - min_freq) / 9, -6)  # Round to nearest 1 MHz

# Generate y-axis tick locations with clean decimal values
clean_decimal_ticks = np.arange(min_freq, max_freq + step, step)
plt.yticks(clean_decimal_ticks)
plt.grid(True)
plt.show()
