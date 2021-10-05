'''
Jäätmed lihatööstuses tihti sisaldavad palju valku. Sellest valgust oleks võimalik toota aminohapped
ja orgaanilised happed, kui hüdrolüüsida jäätmed vees kõrge temperatuuri ja rõhu juures. Kuid katsed
on näidanud, et aminohappe saagis on tavaliselt palju väiksem, kui teoreetiliselt peaks olema.
Näiteks, 1 g kala jäätmetes (kuiv) oli 600 mg valku, aga toodeti ainult 140 mg aminohappeid. Saagis
on madal suuresti sellepärast, et aminohapped ise lagunevad reaktori tingimustel. Suurema saagise
saavutamiseks oleks vaja vähendada aminohapete lagunemist.

Jäätmeid töödeldatakse torureaktoris ja toodetud aminohapped hakkavad lagunema. Torureaktor on
2000-liitrine ja jäätmete kogus algses lahuses on 50 g/l. Kuna valk moodustab umbes 60% nendest
jäätmetest ja umbes 14% kala valgust on glütsiin, siis algne glütsiini kontsentratsioon on 56 mmol/l,
kuigi see on veel valgus. Mis on glütsiini kontsentratsioon väljuvas voos? Reaktor töötab 20 MPa ja
250 kraadi juures. Mahtkulu on 2 l/s.

Hinnangulised kineetilised parameetrid kala jäätmete lagunemise reaktsiooniks:
A = 1*10^12 1/s
Ea = 150 kJ/mol

Kineetilised parameetrid glütsiini lagunemise reaktsiooniks:
A = 3,51*10^13 1/s
Ea = 166 kJ/mol

Viited:
* N. Sato, A. T. Quitain, K. Kang, H. Daimon, and K. Fujie, “Reaction Kinetics of Amino Acid
Decomposition in High-Temperature and High-Pressure Water,” Ind. Eng. Chem. Res., vol. 43, no. 13,
pp. 3217–3222, Jun. 2004, doi: 10.1021/ie020733n.
* K. Kang et al., “Optimization of amino acids production from waste fish entrails by hydrolysis
in sub and supercritical water,” The Canadian Journal of Chemical Engineering, vol. 79, no. 1,
pp. 65–70, 2001, doi: 10.1002/cjce.5450790110.
* B. Mohanty et al., “Amino Acid Compositions of 27 Food Fishes and Their Importance in Clinical
Nutrition,” J Amino Acids, vol. 2014, p. 269797, 2014, doi: 10.1155/2014/269797.
'''
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

RGAS = 8.314462618 # J mol^-1 K^-1

def kiiruskonstant1(T):
    # valgu lagunemise reaktsioon
    A = 1e12
    Ea = 150 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15)) # s^-1
    return k

def kiiruskonstant2(T):
    # aminohappe lagunemise reaktsioon
    A = 3.51e13
    Ea = 166 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15)) # s^-1
    return k

Cjäägid_alg = 0.056 # mol l^-1
V = 2000 # l

m = GEKKO()

T = m.Param(value=250) # Celcius
k1 = m.Param(value=kiiruskonstant1(T.value))
k2 = m.Param(value=kiiruskonstant2(T.value))
vdot = m.Param(value=2) # l s^-1

npts = int(V)
m.time = np.linspace(0, V, npts) # tegelikult siin on maht, mitte aeg. Gekkos saab kasutada ".time" ka mahu tuletise jaoks

Cjäägid = m.Var()
Ca = m.Var()
Cjäägid.value = Cjäägid_alg * np.ones(npts) # alustame algkontsentratsiooniga
Ca.value = np.zeros(npts)

# Meil on mudelis kaks ainet (jäägid ja glütsiin), siis saame koostada kahte moolbilanssi:
m.Equation(Cjäägid.dt() == -k1 / vdot * Cjäägid)
m.Equation(Ca.dt() == k1 / vdot * Cjäägid - k2 / vdot * Ca)

m.options.IMODE = 4 # Gekkos on 6 arvutamise viisi, IMODE=4 tähendab, et lahendame dünaamilist mudelit
m.solve(disp=False)

plt.rcParams.update({'font.size': 13})

plt.figure()
plt.plot(m.time, Cjäägid.value, label="Jäägid")
plt.plot(m.time, Ca.value, label="Glütsiin")
plt.xlabel('Maht (l)')
plt.ylabel('Kontsentratsioon (mol/l)')
plt.legend(frameon=False)
plt.show()
