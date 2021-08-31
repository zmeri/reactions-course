'''
KAT0141 - Reaktsiooniprotsessid

Näide: Reaktsiooni kineetika parameetrite määramine
---------------------------------------------------
Reaktsioon triphenüül metüül kloriidi (A) ja metanooli (B) vahel viiakse läbi benseeni ja
püridiini lahuses 25 kraadi juures. Püridiin reageerib HCl'iga ja sadeneb välja
kui püridiin vesinikkloriid, ja selle tõttu reaktsioon on ühesuunaline.

Perioodilises reaktoris mõõdeti A kontsentratsioon ajas. Andmed antakse failis
harjutustund2_1_andmed.csv. Metanooli algkontsentratsioon oli 0.5 mol dm^-3.

Osa 1: Leia reaktsiooni järk triphenüül metüül kloriidi suhtes.
Osa 2: Teistes katsetes leiti, et reaktsiooni järk metanooli suhtes oli 1. Arvuta
    reaktsiooni kiirusekonstanti.

Viide: Näide 5-1 raamatus Fogler, H.S. Elements of Chemical Reaction Engineering, 4th Ed. Pearson Education, Inc, 2006.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import minimize

# Andmete import ---------------------------------------------------------------
filename = 'kineetika_määramise_andmed.csv'
col_names = ['aeg', 'kont']
df = pd.read_csv(filename, skiprows=2, delimiter=';', names=col_names)

print(df)

plt.figure()
plt.scatter(df['aeg'], df['kont'])
plt.title('Toorandmed')
plt.xlabel('Aeg (min)')
plt.ylabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.show()

"""
Reaktsiooni kiiruseseadus:
-r_A = k * C_A**a * C_B**b
aga kuna C_B alguses on palju kõrgem kui C_A võime eeldata, et C_B ei muutu
-r_A = k * C_A**a * C_B0**b
-r_A = (k * C_B0**b) * C_A**a
-r_A = k' * C_A**a
"""

# Tuletise võtmine ------------------------------------------------------------
df = df.assign(dc_dt = np.gradient(df['kont'], df['aeg'], edge_order=2))
print(df)

plt.figure()
plt.scatter(df['kont'], df['dc_dt'])
plt.title('Tuletis (dC_dt)')
plt.xlabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.ylabel('Kiirus (mol dm$^{-3}$ min$^{-1}$)')
plt.show()

# Lineaar regressioon -----------------------------------------------------------------
tulemus = linregress(np.log(df['kont'].to_numpy()), np.log(-df['dc_dt'].to_numpy()))
print('regressiooni tulemus:\n', tulemus)

kont_grid = np.linspace(df['kont'].min(), df['kont'].max(), 50)
dc_dt_calc = np.exp(tulemus.slope * np.log(kont_grid) + tulemus.intercept)

print('\nreaktsiooni järk=', tulemus.slope)

plt.figure()
plt.scatter(df['kont'], -df['dc_dt'])
plt.plot(kont_grid, dc_dt_calc)
plt.title('Lineaar regressioon')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log(C$_A$)')
plt.ylabel('log(-dc/dt)')
plt.show()

# Regressioon --------------------------------------------------------------------
def kiiruseseadus_reg(par, c_a, r):
    r_calc = par[0] * c_a**par[1]
    return np.sum(np.log(r_calc/r)**2)
    # return np.sum((r_calc - r)**2)

def kiiruseseadus(par, c_a):
    return par[0] * c_a**par[1]

guess = np.asarray([0.13, 2])
bnds = ((1e-5, 0.5), (0,3))
tulemus = minimize(kiiruseseadus_reg, guess, bounds=bnds, args=(df['kont'].to_numpy(), -df['dc_dt'].to_numpy()))
print('regressiooni tulemus:\n', tulemus)

kont_grid = np.linspace(df['kont'].min(), df['kont'].max(), 50)
dc_dt_calc = kiiruseseadus(tulemus.x, kont_grid)

print('\nreaktsiooni järk=', tulemus.x[1])

plt.figure()
plt.scatter(df['kont'], -df['dc_dt'])
plt.plot(kont_grid, dc_dt_calc)
plt.title('Regressioon lahendajaga')
plt.xlabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.ylabel('-dc/dt')
plt.show()

# Arvuta kiirusekonstanti ----------------------------------------------------
c_b0 = 0.5 # kontsentratsioon B alguses (mol dm^-3)
k = tulemus.x[0] / c_b0
print('\nk=', k, '(dm^3 mol^-1)^2 min^-1')

# graafiline meetod ------------------------------------------------------------
# null järk
plt.figure()
plt.scatter(df['aeg'], df['kont'])
plt.xlabel('Aeg (min)')
plt.ylabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.title('0. järk')
plt.show()

# esimene järk
plt.figure()
plt.scatter(df['aeg'], np.log(df['kont']))
plt.xlabel('Aeg (min)')
plt.ylabel('log(C$_A$)')
plt.title('1. järk')
plt.show()


# teine järk
plt.figure()
plt.scatter(df['aeg'], 1/df['kont'])
plt.xlabel('Aeg (min)')
plt.ylabel('1/C$_A$')
plt.title('2. järk')
plt.show()
