import numpy as np

epsilon_s = (8.85e-14)*11.7
q = 1.602e-19

epsilon_ox = 3.9*(8.85e-14)
const = (1/0.3)**2-1

C_ox = np.sqrt(const/3*q*epsilon_s*2.15*10e17)
print(C_ox)

tox = epsilon_ox/C_ox
print(tox)

C_ox = 3.467e-7
xdo = epsilon_s/C_ox*(np.sqrt(1+(2*C_ox**2)/(epsilon_s*q*2.1410e17))-1)
print(xdo)

phi_gb0 = 1/2*q*2.14e17/epsilon_s*xdo**2
print(phi_gb0)

def xd(v_gb):
    return epsilon_s/C_ox*(np.sqrt(1+(2*C_ox**2*(1+v_gb))/(epsilon_s*q*2.1410e17))-1)
def surface_potential(v_gb):
    return q*2.14e17/(2*epsilon_s)*xd(v_gb)**2
print(xd(0.5))
print(surface_potential(0.5))
