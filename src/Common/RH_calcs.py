import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scc

file = "Jouan_initial_condition.txt"
matrix = np.loadtxt(file, skiprows=1)

T = matrix[:, 0]   # Temp (K)
Td = matrix[:, 1]   # Dewpoint (K)
Qv = matrix[:, 2]   # Raw Qv from file (kg/kg)
Qsat = matrix[:, 3]   # Qsat from file (kg/kg)
Pressure = matrix[:, 4] # Pascals
alt = matrix[:, 9]   # Altitude (meters)

T_fahr = scc.convert_temperature(T, "K", "F")
Td_fahr = scc.convert_temperature(Td, "K", "F")

T_c = scc.convert_temperature(T, "K", "C")
Tdc = scc.convert_temperature(Td, "K", "C")

es = 6.1078 * np.exp((17.269*T_c)/(237.3 + T_c))
ea = 6.1078 * np.exp((17.269*Tdc)/(237.3+Tdc))

plt.plot(Qv/Qsat * 100, alt, label="Qsat depend (Qv/Qsat * 100)")
plt.plot(100 - 5*(T - Td), alt, label="Temp depend (100 - 5*(T-Td))")
plt.plot(100*((112-(0.1*T_fahr) + Td_fahr)/(112 + (0.9*T_fahr))**8), label="found online")
plt.plot(ea/es * 100, alt, label="Magnus tetens")
plt.ylabel("altitude")
plt.xlabel("RH")
plt.title("Relative humidity?")
plt.legend()
plt.show()