# impordime vajalikud paketid
import matplotlib.pyplot as plt
import numpy as np
from gekko import GEKKO

# defineerime lähteandmeid
m = GEKKO()
CA0 = m.Param(value=0.5) # mol/l
CB0 = m.Param(value=0.5) # mol/l
k = m.Param(value=0.008) # l mol^-1 s^-1

t_lõpp = 600 # s
npts = 200
m.time = np.linspace(0, t_lõpp, npts)

# määrame muutujaid
XA = m.Var(0)

# defineerime võrrandeid
m.Equation(XA.dt() == k * (1 - XA) * (CB0 - CA0 * XA))

# kutsume lahendaja
m.options.IMODE = 4
m.solve(disp=False)

# väljastame tulemust
plt.figure()
plt.plot(m.time, XA.value)
plt.xlabel('Aeg (s)')
plt.ylabel('Konversiooniaste')
plt.show()
