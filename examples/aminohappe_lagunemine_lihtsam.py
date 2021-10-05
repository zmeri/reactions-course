'''
Jäätmed lihatööstuses tihti sisaldavad palju valku. Sellest valgust oleks võimalik toota aminohapped
ja orgaanilised happed, kui hüdrolüüsida jäätmed vees kõrge temperatuuri ja rõhu juures. Kuid katsed
on näidanud, et aminohappe saagis on tavaliselt palju väiksem, kui teoreetiliselt peaks olema. Näiteks,
1 g kalajäätmetes (kuiv) oli 600 mg valku, aga toodeti ainult 140 mg aminohappeid. Saagis on madal
suuresti sellepärast, et aminohapped ise lagunevad reaktori tingimustel. Suurema saagise saavutamiseks
oleks vaja vähendada aminohapete lagunemist.

Jäätmeid töödeldatakse torureaktoris ja toodetud aminohapped hakkavad lagunema. Selleks, et lihtsustada ülesannet keskendume ainult aminohapete lagunemisele. Eeldame, et torureaktor on 120-liitrine ja et aminohappeid juurde ei teki. Alguses ühe aminohappe (glütsiini) kontsentratsioon on 10 mmol/l. Mis osa glütsiinist on lagunenud selleks ajaks, kui lahus on väljunud reaktorist? Reaktor töötab 20 MPa ja 250
kraadi juures ja mahtkulu on 2 l/s.

Kineetilised parameetrid glütsiini lagunemiseks:
A = 3,51*10^13 1/s
Ea = 166 kJ/mol

Viide: N. Sato, A. T. Quitain, K. Kang, H. Daimon, and K. Fujie, “Reaction Kinetics of Amino Acid
Decomposition in High-Temperature and High-Pressure Water,” Ind. Eng. Chem. Res., vol. 43, no. 13,
pp. 3217–3222, Jun. 2004, doi: 10.1021/ie020733n.
'''
import numpy as np

RGAS = 8.314462618 # J mol^-1 K^-1

def kiiruskonstant(T):
    A = 3.51e13
    Ea = 166 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15)) # s^-1
    return k

Ca_alg = 0.01 # mol l^-1
V = 120 # l
T = 250 # Celcius
vdot = 2 # l s^-1
ndot_a_alg = vdot * Ca_alg

k = kiiruskonstant(T)

'''
Torureaktori põhivõrrandis on vaja lahendada integraali konversiooniastme suhtes. Kui teha seda
saame järgmine valem:
    V = -ndot_a_alg / Ca_alg / k * ln(1 - Xa)
Lahendades Xa jaoks valem muutub selliseks:
'''
Xa = 1 - np.exp(-V * Ca_alg * k / ndot_a_alg)
print('Aminohappe konversiooniaste= {:.2f}%'.format(Xa * 100))
