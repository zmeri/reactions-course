'''
Epikloorhüdriin on oluline lähteaine epoksiidvaikude, elastomeeride ja sünteetilise glütseriini tootmiseks. Seda praegu toodetakse enamasti propüleenist ja kloorgaasist, aga on võimalik ka toota glütseroolist. Tootmine glütseroolist on võibolla eelistatud kuna ei ole vaja kasutada kloorgaasi. Esimeses sammus toimub reaktsioon kloorhappe ja glütserooli vahel katalüsaatori juuresolekul selleks, et toota 2,3-dikloro-1-propanooli (DKP). Seejärel viiakse läbi reaktsiooni DKP ja naatrium hüdroksiidi vahel epikloorhüdriini sünteesimiseks.

C3H6OCl2 + NaOH -> C3H5OCl + H2O + NaCl

Tahetakse läbi viia selle reaktsiooni 1-liitrilises perioodilises reaktoris. Kui kaua aega segu peaks olema reaktoris selleks, et saavutada 90% konversiooniastet?

DKP algkontsentratsioon on 0,08 mol/l ja naatrium hüdroksiidi algkontsentratsioon on 0,12 mol/l. Reaktorit hoitakse 60 kraadi juures. DKP suhtes on teist järku reaktsioon ja NaOH suhtes esimest järku. Reaktsiooni kineetilised parameetrid:
A = 1,61*10^25 l^3 mol^-2 s^-1
Ea = 150 kJ/mol

Viide: J. S. Zhang, Y. C. Lu, Q. R. Jin, K. Wang, and G. S. Luo, “Determination of kinetic parameters of dehydrochlorination of dichloropropanol in a microreactor,” Chemical engineering journal, vol. 203, pp. 142–147, 2012.
'''
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

RGAS = 8.314462618 # J mol^-1 K^-1

def kiiruskonstant(T):
    A = 1.61e25 # l^3 mol^-2 s^-1
    Ea = 150 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15))
    return k # l^3 mol^-2 s^-1

m = GEKKO()

C_dkp_alg = m.Param(value=0.08) # mol l^-1
C_oh_alg = m.Param(value=0.12) # mol l^-1
V = m.Param(value=1) # l
T = m.Param(value=60) # Celcius
k = m.Param(value=kiiruskonstant(T.value))

t_lõpp = 200 # s
npts = 1000
m.time = np.linspace(0, t_lõpp, npts)

X_dkp = m.Var() # konversiooniaste
X_dkp.value = np.zeros(npts)

# paneme kirja põhivõrrandit diferentsiaalkujul
m.Equation(X_dkp.dt() == k * C_dkp_alg * (1 - X_dkp)**2 * (C_oh_alg - C_dkp_alg * X_dkp))

m.options.IMODE = 4
m.solve(disp=False)

idx = np.amin(np.where(np.asarray(X_dkp.value) > 0.9)[0])
print('Vajalik aeg= {:.2f} s'.format(m.time[idx]))

plt.figure()
plt.plot(m.time, X_dkp)
plt.xlabel('Aeg (s)')
plt.ylabel('Konversiooniaste')
plt.show()
