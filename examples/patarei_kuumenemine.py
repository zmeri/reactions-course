'''
Kui akud on ebastabiilses seisundis, ühendid akus võivad hakata reageerima, põhjustades tulekahju. Selline olukord võib näiteks tekkida, kui aku temperatuur on liiga kõrge või kui laaditakse aku liiga palju. On võimalik teha turvalisemaid akusid, kui on mudelid, mis kirjeldavad neid protsesse akus, mis võivad viia tulekahjuni.

Oletame, et liitium-ioon aku pannakse keskkonda, mis on 240 kraadi juures. Kui kaua läheks enne, kui reaktsioonid aku sees põhjustavad järsku temperatuuri tõusu (ja tõenäoliselt tulekahju)? Vajalikud parameetrid antakse Hewsoni artiklis.

Viited
------
* J. C. Hewson, “Understanding the limits of thermal runaway in lithium-ion battery systems,” Sandia National Lab. (SNL-NM), Albuquerque, NM (United States), SAND2016-6054C, Jun. 2016. Accessed: Oct. 21, 2021. [Online]. Available: https://www.osti.gov/biblio/1420927
'''
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

m = GEKKO()

# Mudelis on neli reaktsiooni:
# 1) Metastabiilse tahke-elektrolüüdi liidese (ingl solid-electrolyte interface, SEI) lagunemine
# 2) Anoodi-elektrolüüdi reaktsioon
# 3) Elektrolüüdi lagunemine
# 4) CoO2-elektrolüüdi reaktsioon

# parameetrid artiklist
Yalg_väärtused = np.asarray([0.15, 0.45, 1.0, 0.96]) # Y on sisuliselt normaliseeritud kontsentratsioon, ehk näitab, kui palju antud ainest on järgi jäänud
A_väärtused = np.asarray([1.67e15, 1.67e6, 5.1e25, 6.67e11]) # s^-1, sagedustegur
Ea_väärtused = np.asarray([135, 77.2, 274, 122]) * 1000 # J mol^-1, aktivatsioonienergia
Hrxn_väärtused = np.asarray([257.0, 1714, 155, 314]) * 1000 # J kg^-1, reaktsioonientalpia
rho_väärtused = np.asarray([100, 610.0, 1221, 407]) # kg m^-3, aine mass mahu ühiku kohta

ncomp = 4 # reaktsioonide arv
Yalg = m.Array(m.Param, (ncomp))
A = m.Array(m.Param, (ncomp))
Ea = m.Array(m.Param, (ncomp))
Hrxn = m.Array(m.Param, (ncomp))
rho = m.Array(m.Param, (ncomp))
for i in range(ncomp):
    Yalg[i].value = Yalg_väärtused[i]
    A[i].value = A_väärtused[i]
    Ea[i].value = Ea_väärtused[i]
    Hrxn[i].value = Hrxn_väärtused[i]
    rho[i].value = rho_väärtused[i]

RGAS = m.Const(value=8.314462618) # J mol^-1 K^-1
rhoCp = m.Param(value=2.79e6) # J m^-3
pind_maht_suhe = m.Param(48.4) # m^-1
epsilon = m.Param(0.8)
sigma = m.Param(5.67e-8) # W m^-2 K^-4
lm = 0.034 # W m^-1 K^-1
d = 0.018 # m
Nu = 4.8
hc = Nu * lm / d # W m^-2 K^-1
hbl = 7.17 # W m^-2 K^-1, soojusläbikandetegur piirkihi jaoks
heff = m.Param(value = 1 / (1 / hc + 1 / hbl)) # soojusläbikandetegur
Tvälis = m.Param(value=240) # K, ümbritseva keskkonna temperatuur

npts = 1000
t_lõpp = 3400
m.time = np.linspace(0, t_lõpp, npts) # s

T = m.Var(25)
Y0 = m.Var(Yalg[0].value)
Y1 = m.Var(Yalg[1])
Y2 = m.Var(Yalg[2])
Y3 = m.Var(Yalg[3])

z = m.Intermediate(0.033 + (Yalg[1] - Y1) + (Yalg[0] - Y0))
k = [None] * ncomp
for i in range(ncomp):
    k[i] = m.Intermediate(A[i] * m.exp(-Ea[i] / RGAS / (T + 273.15)))

Qrxn0 = m.Intermediate(k[0] * Y0 * rho[0] * Hrxn[0])
Qrxn1 = m.Intermediate(k[1] * Y1 * m.exp(-z/0.033) * rho[1] * Hrxn[1])
Qrxn2 = m.Intermediate(k[2] * Y2 * rho[2] * Hrxn[2])
Qrxn3 = m.Intermediate(k[3] * Y3 * (1 - Y3) * rho[3] * Hrxn[3])
Qkonv = m.Intermediate(heff * pind_maht_suhe * (Tvälis - T))
Qkiirg = m.Intermediate(epsilon * sigma * pind_maht_suhe * ((Tvälis + 273.15)**4 - (T + 273.15)**4))
m.Equation(T.dt() == (Qrxn0 + Qrxn1 + Qrxn2 + Qrxn3 + Qkonv + Qkiirg) / rhoCp)
m.Equation(Y0.dt() == -k[0] * Y0)
m.Equation(Y1.dt() == -k[1] * Y1 * m.exp(-z/0.033))
m.Equation(Y2.dt() == -k[2] * Y2)
m.Equation(Y3.dt() == -k[3] * Y3 * (1 - Y3))

m.options.IMODE = 4
m.solve(disp=False)

plt.figure()
plt.plot(m.time, T.value)
plt.xlabel('Aeg (s)')
plt.ylabel("Temperatuur (°C)")
plt.show()
