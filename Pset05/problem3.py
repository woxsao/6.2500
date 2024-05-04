import numpy as np

epsilon_s = (8.85e-14)*11.7
q = 1.602e-19

epsilon_ox = 3.9*(8.85e-14)

t_ox = 2.5e-7
Na = 1e18
Nd = 1e19
L = 110e-7
ni = 1e10
mu_n = 775
W = 2e-4
eta = 1
Vdd = 1


psi_p = 0.025*np.log(Na/ni)
print("psi_p", psi_p)
V_fb = -1*(0.56+(0.06*np.log10(Na/ni)))
print("V_fb", V_fb)
C_ox = epsilon_ox/t_ox

V_t = V_fb + 2*psi_p + 1/C_ox*np.sqrt(2*epsilon_s*q*Na*2*psi_p)
print("V_t", V_t)

K_on = mu_n*C_ox/(2*L)
print("K_on", K_on)
I_on = W*K_on*(Vdd-V_t)**2
print("I_on", I_on)

Dn = 0.025*mu_n
Ld = np.sqrt(2*epsilon_s*0.025/(q*Na))
Q_no = q*Na*np.sqrt(np.pi)/2*Ld
K_off = Q_no/L*Dn

print("K_off", K_off)

I_off = W*K_off*np.exp(-V_t/(eta*0.025))
print("I_off", I_off)


C_in = epsilon_ox*L*W*2/t_ox
Cw = 9e-15
C_out = 1/4*C_in
print("Cout", C_out)
print("Cin", C_in)
alpha = 0.1
Ed = (C_out+Cw+2*C_in)*Vdd**2*alpha
print("Ed", Ed)
ngate_per_stage = 2**18-1
nstages = 5


FO = 2
LD = 18
Rw = 1200
tclk = LD*1/2*np.log(2)*((Vdd/I_on)*(C_out + Cw+2*C_in)+Rw*(Cw/2+FO*C_in))


Ed_per_clock = ngate_per_stage*nstages*Ed
print("Ed_per_clock", Ed_per_clock)
E_leak = Vdd*I_off*tclk
print("Eleak", E_leak)
E_leak_per_clock = ngate_per_stage*nstages*E_leak
print("Eleak per clock", E_leak_per_clock)

print("SUMMARY---------------------------")
print("FOUND TOTAL ENERGY PER CLOCK", E_leak_per_clock + Ed_per_clock)
print("REQUIRED ENERGY PER CLOCK", 3e-9)
print("FOUND TCLK", tclk)
print("REQUIRED tclk", 1/1.7e9)