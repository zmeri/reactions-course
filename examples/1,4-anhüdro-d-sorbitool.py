'''
Praegu keemiatööstus kasutab suures osas naftat toorainena. Selleks, et muuta majandust jätkusuutlikumaks, otsitakse viise kemikaalide tootmiseks taastuvatest allikatest. Biomass võib olla üks nendest taastuvatest allikatest. Sorbitool on üks kemikaal, mida on võimalik tselluloosist toota. Sorbitoolist võib omakorda toota veel mitu produkti. Näiteks, on leitud, et on võimalik toota 1,4-anhüdro-d-sorbitooli, kui läbi viia reaktsiooni vees kõrgel temperatuuril. Kuid toimuvad samaaegselt veel mitu reaktsiooni, mis vähendavad 1,4-anhüdro-d-sorbitooli saagist.

Kui on 10 liitrine segureaktor, mis peaks olema siseneva voo mahtkulu selleks, et 1,4-anhüdro-d-sorbitooli saagis oleks 40%?

Reaktor töötab temperatuuril 275 °C ja sorbitooli kontsentratsioon sisenevas voos on 1,5 mol/l. Viidetud artiklis antakse reaktsioonide Arrheniuse võrrandi parameetrid.

Viide
-----
* Yamaguchi A, Hiyoshi N, Sato O, Shirai M. Sorbitol dehydration in high temperature liquid water. Green Chemistry. 2011;13(4):873-81.
'''
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

# A = d-sorbitool
# B = 1,4-anhüdro-d-sorbitool
# C = 2,5-anhüdro-d-sorbitool
# D = isosorbiid

nrxn = 4 # reaktsioonide arv
CA0 = 1.5 # mol l^-1
saagis = 0.4
CB = saagis * CA0
print('CB = ', CB)

m = GEKKO()

Ea_numpy = np.asarray([127000., 166000., 195000, 136000]) # J mol^-1
A_numpy = np.asarray([6.2358e+11, 4.8870e+14, 4.7301e+17, 3.1527e+11]) # h^-1

Ea = m.Array(m.Param, (nrxn))
A = m.Array(m.Param, (nrxn))
for i in range(nrxn):
    Ea[i].value = Ea_numpy[i]
    A[i].value = A_numpy[i]

T = m.Param(value=275) # degC reaktori temperatuur
CA0 = m.Param(value=CA0) # mol l^-1
CB = m.Param(value=CB) # mol l^-1
V = m.Param(value=10) # l
RGAS = m.Const(value=8.314) # J mol^-1 K^-1

vdot = m.Var()
CA = m.Var()
CC = m.Var()
CD = m.Var()

# arvuta reaktsioonide kiirusekonstanti
k = [None] * nrxn
for i in range(nrxn):
    k[i] = m.Intermediate(A[i] * m.exp(-Ea[i] / RGAS / (T + 273.15)))

# mudeli valemid
m.Equation(0 == vdot * CA0 - vdot * CA + (-k[0] * CA - k[1] * CA) * V)
m.Equation(0 == -vdot * CB + (k[0] * CA - k[2] * CB) * V)
m.Equation(0 == -vdot * CC + k[1] * CA * V)
m.Equation(0 == -vdot * CD + (k[1] * CB - k[3] * CD) * V)

# lahendada ja printida vastuseid
m.solve(disp=False)

print('Mahtkulu= {} l/h'.format(vdot.value[0]))
print('CA= {} mol/l'.format(CA.value[0]))
print('CC= {} mol/l'.format(CC.value[0]))
print('CD= {} mol/l'.format(CD.value[0]))
