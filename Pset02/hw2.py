import numpy as np
from sympy.solvers import solve
from sympy import Symbol

#Problem 2a
Nd = 10**16
ni = 10**10
E_g = 1.12
k = 8.62*10**-5
T = 300
A_si = 4.7*10**15
q = 1.60217663 * 10**19

efei = k*T*np.log(Nd/ni)
print("ef-ei: ", efei)

#problem 2bv
Na = 10**16
Nd = 10**16
print("phi_b", 2*efei)


#problem 3b
ef_prime = k*300*np.log(100)
print("new ef with 100x conductance", ef_prime)
#problem 3c
T_new = 363.15
ni_new = A_si*T_new**(3/2)*np.exp(-E_g/(2*k*T_new))
print("ni_new at 90C", ni_new)
V = k*T_new*np.log(100*10**16/ni_new)
print("V at 90C:", V)

T_new = 77
ni_new = A_si*T_new**(3/2)*np.exp(-E_g/(2*k*T_new))
print("ni_new at 77 K", ni_new)
V = k*T_new*np.log(100*10**16/ni_new)
print("V at 77K:", V)

T_new = 300
ni_new = A_si*T_new**(3/2)*np.exp(-E_g/(2*k*T_new))
print("ni_new at {T_new} K", ni_new)
V = k*T_new*np.log(100*10**16/ni_new)
print("V at {T_new}K:", V)

#problem 3g
phi_b = 2*k*T_new*np.log(10**6)*1000
print("phi_b", phi_b)
Na = np.exp((phi_b+100)/25)*10**4
print("Na", Na)

print("K*T at 90C", k*363.15*np.log(10))
print("K*T at 25C", k*300*np.log(10))
print("K*T at 77K", k*77*np.log(10))

def Map(self, keyvalue, value):
        results = []
        i = 0
        n = len(value)
        while i < n:
            # skip non-ascii letters in C/C++ style a la MapReduce paper:
            while i < n and value[i] not in string.ascii_letters:
                i += 1
            start = i
            while i < n and value[i] in string.ascii_letters:
                i += 1
            w = value[start:i]
            if start < i and w.istitle():
                results.append ((w.lower(), 1))
        return results
