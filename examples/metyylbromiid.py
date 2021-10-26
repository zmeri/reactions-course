'''
Metüülbromiidi  saamiseks  viiakse  läbi  pöördumatu  vedelfaasiline
lihtreaktsioon poolperioodilises reaktoris:
    CNBr + CH3NH2 ---> CH3Br + NCNH2
    ehk  A  +  B ----> C  +  D
Alguses on juba reaktoris 5 liitrit lahust, mis sisaldab 0,05 mol/l ainet
A. Siis hakatakse lisama lahust, mis sisaldab aine B kontsentratsioonis
0,025 mol/l. Leida  broomtsüaniidi  konversiooniastme sõltuvus ajast.
Samuti leida broomtsüaniidi, metüülamiini  ja  metüülbromiidi
kontsentratsiooni sõltuvus ajast.

Teist lahust lisatakse mahtkuluga 0,05 l s^-1. Reaktsiooni kiiruskonstant
antud tingimustel on 2,2 l s^-1 mol^-1.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

m = GEKKO()

CA0 = m.Param(value=0.05) # mol/L
CB0 = m.Param(value=0.025) # mol/L
V0 = m.Param(value=5) # L
k = m.Param(value=2.2) # L/s/mol
v_in = 0.05 # L/s
nA0 = m.Param(value=CA0.value * V0.value)

npts = 1000
t_lõpp = 400
m.time = np.linspace(0, t_lõpp, npts) # s

V = m.Param(value=V0.value + v_in * np.asarray(m.time))
FB0 = m.Param(value=CB0.value * v_in * np.asarray(m.time))

XA = m.Var(0)

m.Equation(XA.dt() == k * (1 - XA) * (FB0 - nA0 * XA) / V)

m.options.IMODE = 4
m.solve(disp=False)

CB = (FB0.value - nA0.value * np.asarray(XA.value)) / np.asarray(V.value)
CC = (nA0.value * np.asarray(XA.value)) / np.asarray(V.value)
CA = nA0.value * (1 - np.asarray(XA.value)) / np.asarray(V.value)

plt.figure()
plt.plot(m.time, XA.value)
plt.xlabel('Aeg (s)')
plt.ylabel('X$_A$')
plt.savefig('../media/metyylbromiid_XA.png')
plt.show()

plt.figure()
plt.plot(m.time, CA, label='CA')
plt.plot(m.time, CB, label='CB')
plt.plot(m.time, CC, label='CC')
plt.xlabel('Aeg (s)')
plt.ylabel('Kontsentratsioon (mol L$^{-1}$)')
plt.legend(frameon=False)
plt.savefig('../media/metyylbromiid_kontsentratsioonid.png')
plt.show()
