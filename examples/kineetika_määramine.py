'''
KAT0141 - Reaktsiooniprotsessid

Näide: Reaktsiooni kineetika parameetrite määramine
---------------------------------------------------
Reaktsioon triphenüül metüül kloriidi (A) ja metanooli (B) vahel viiakse läbi benseeni ja
püridiini lahuses 25 kraadi juures. Reaktsioon on ühesuunaline kuna püridiin reageerib HCl'iga
ja sadeneb välja kui püridiin vesinikkloriid.

Perioodilises reaktoris mõõdeti A kontsentratsioon ajas. Andmed antakse failis
kineetika_määramise_andmed.csv. Metanooli algkontsentratsioon oli 0,5 mol dm^-3.

Osa 1: Leia reaktsiooni järk triphenüül metüül kloriidi suhtes.
Osa 2: Teistes katsetes leiti, et reaktsiooni järk metanooli suhtes oli 1. Arvuta
    reaktsiooni kiiruskonstanti.

Viide: Näide 5-1 raamatus Fogler, H.S. Elements of Chemical Reaction Engineering, 4th Ed. Pearson Education, Inc, 2006.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import minimize

# Andmete import ---------------------------------------------------------------
filename = 'kineetika_määramise_andmed.csv'
col_names = ['aeg', 'kontA']
df = pd.read_csv(filename, skiprows=2, delimiter='\t', names=col_names)

print(df)

plt.figure()
plt.scatter(df['aeg'], df['kontA'])
plt.title('Toorandmed')
plt.xlabel('Aeg (min)')
plt.ylabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.show()

"""
Reaktsiooni kiirusevalem:
-r_A = k * C_A**a * C_B**b
aga kuna C_B alguses on palju kõrgem kui C_A võime eeldata, et C_B ei muutu
-r_A = k * C_A**a * C_B0**b
-r_A = (k * C_B0**b) * C_A**a
-r_A = k' * C_A**a
"""

# Tuletise võtmine ------------------------------------------------------------
df['dc_dt'] = np.gradient(df['kontA'], df['aeg'], edge_order=2)
print(df)

plt.figure()
plt.scatter(df['kontA'], df['dc_dt'])
plt.title('Tuletis (dC_dt)')
plt.xlabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.ylabel('Kiirus (mol dm$^{-3}$ min$^{-1}$)')
plt.show()

# Lineaar regressioon -----------------------------------------------------------------
tulemus = linregress(np.log(df['kontA'].to_numpy()), np.log(-df['dc_dt'].to_numpy()))
print('regressiooni tulemus:\n', tulemus)

kont_grid = np.linspace(df['kontA'].min(), df['kontA'].max(), 50)
dc_dt_calc = np.exp(tulemus.slope * np.log(kont_grid) + tulemus.intercept)

print('\nreaktsiooni järk=', tulemus.slope)

plt.figure()
plt.scatter(df['kontA'], -df['dc_dt'])
plt.plot(kont_grid, dc_dt_calc)
plt.title('Lineaar regressioon')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log(C$_A$)')
plt.ylabel('log(-dc/dt)')
plt.show()

# Regressioon --------------------------------------------------------------------
def kiirusevalem_reg(par, c_a, r):
    r_calc = par[0] * c_a**par[1]
    return np.sum(np.log(r_calc/r)**2)

def kiirusevalem(par, c_a):
    return par[0] * c_a**par[1]

guess = np.asarray([0.13, 2])
bnds = ((1e-5, 0.5), (0,3))
tulemus = minimize(kiirusevalem_reg, guess, bounds=bnds, args=(df['kontA'].to_numpy(), -df['dc_dt'].to_numpy()))
print('regressiooni tulemus:\n', tulemus)

kont_grid = np.linspace(df['kontA'].min(), df['kontA'].max(), 50)
dc_dt_calc = kiirusevalem(tulemus.x, kont_grid)

print('\nreaktsiooni järk=', tulemus.x[1])

plt.figure()
plt.scatter(df['kontA'], -df['dc_dt'])
plt.plot(kont_grid, dc_dt_calc)
plt.title('Regressioon lahendajaga')
plt.xlabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.ylabel('-dc/dt')
plt.show()

# Arvuta kiiruskonstanti ----------------------------------------------------
c_b0 = 0.5 # B kontsentratsioon alguses (mol dm^-3)
k = tulemus.x[0] / c_b0
print('\nk=', k, '(dm^3 mol^-1)^2 min^-1')

# graafiline meetod ------------------------------------------------------------
# null järk
plt.figure()
plt.scatter(df['aeg'], df['kontA'])
plt.xlabel('Aeg (min)')
plt.ylabel('Kontsentratsioon (mol dm$^{-3}$) x 1000')
plt.title('0. järk')
plt.show()

# esimene järk
plt.figure()
plt.scatter(df['aeg'], np.log(df['kontA']))
plt.xlabel('Aeg (min)')
plt.ylabel('log(C$_A$)')
plt.title('1. järk')
plt.show()


# teine järk
plt.figure()
plt.scatter(df['aeg'], 1/df['kontA'])
plt.xlabel('Aeg (min)')
plt.ylabel('1/C$_A$')
plt.title('2. järk')
plt.show()
