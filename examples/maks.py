'''
Kui alkoholi juuakse, siis keha reageerib ja proovib seda etanooli muundata selleks, et vähendada etanooli kontsentratsiooni kehas. Enamasti, etanooli muundamine toimub maksas. Etanool maos ja soolastikus imetakse verre, ja kui veri jõuab maksani ensüümid hakkavad seda muundama.

Etanooli muundatakse maksas järgnevate reaktsioonide kaudu:
    CH3CH2OH + NAD+ <-> CH3CHO + NADH + H+
    CH3CHO + NAD+ -> CH3COOH + NADH + H+
Ensüümid kiirendavad neid reaktsioone. Esimeses reaktsioonis alkoholdehüdrogenaas muudab etanooli atsetaldehüüdiks ja teises reaktsioonis aldehüüddehüdrogenaase abil atsetaldehüüd muutub atsetaadiks.

Tee mudelit sellest, kuidas etanooli kontsentratsioon muutub peale seda, kui juuakse alkohoolset jooki. Kui juuakse 400 ml veini (etanooli sisaldus 12 vol%), kas etanooli kontsentratsioon veres ületaks Eesti alkoholi piirmäära (0,2 promilli, ehk 0,2 mg/g)? Kui ületab, kui kaua läheks enne, kui etanooli kontsentratsioon veres on jälle madalam, kui piirmäär?
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

m = GEKKO()

# lähteandmed
V_jook = 400 # ml, joogi maht
mahtosa_etanool = 0.12
etanooli_tihedus = 0.78954976 # g/ml
m_etanool_sisse = V_jook * mahtosa_etanool * etanooli_tihedus

maksu_mass = 1.5 # kg, Allikad: C. Eipel, K. Abshagen, and B. Vollmar, “Regulation of hepatic blood flow: The hepatic arterial buffer response revisited,” World Journal of Gastroenterology : WJG, vol. 16, no. 48, p. 6046, Dec. 2010, doi: 10.3748/wjg.v16.i48.6046. AND https://hypertextbook.com/facts/2004/MaryPennisi.shtml
mm = 46.07 # etanooli molaar mass (g/mol)
vere_tihedus = 1.0506 # g/ml, Allikas: R. J. Trudnowski and R. C. Rico, “Specific Gravity of Blood and Plasma at 4 and 37 °C,” Clinical Chemistry, vol. 20, no. 5, pp. 615–616, May 1974, doi: 10.1093/clinchem/20.5.615.
maksu_tihedus = 1.051 # g/ml, Allikas: B. A. Overmoyer, C. E. McLaren, and G. M. Brittenham, “Uniformity of liver density and nonheme (storage) iron distribution,” Arch Pathol Lab Med, vol. 111, no. 6, pp. 549–554, Jun. 1987.
Vmax = m.Param(value=10) # min^-1, Allikas: F. W. Wagner, A. R. Burger, and B. L. Vallee, “Kinetic properties of human liver alcohol dehydrogenase: oxidation of alcohols by class I isoenzymes,” Biochemistry, vol. 22, no. 8, pp. 1857–1863, Apr. 1983, doi: 10.1021/bi00277a018.
Km = m.Param(value=1.2/1000) # M, Allikas: F. W. Wagner, A. R. Burger, and B. L. Vallee, “Kinetic properties of human liver alcohol dehydrogenase: oxidation of alcohols by class I isoenzymes,” Biochemistry, vol. 22, no. 8, pp. 1857–1863, Apr. 1983, doi: 10.1021/bi00277a018.
V_maks = m.Param(value=maksu_mass / maksu_tihedus) # l
V_veri = m.Param(value=5) # l, Allikas: J. Feldschuh and Y. Enson, “Prediction of the normal blood volume. Relation of blood volume to body habitus.,” Circulation, vol. 56, no. 4, pp. 605–612, Oct. 1977, doi: 10.1161/01.CIR.56.4.605.
V_magu = m.Param(value=0.5 + V_jook/1000) # l, mao maht
V_sool = m.Param(value=0.06) # l, soolestiku maht
vdot_magu_sool = m.Param(value=0.01) # l/min, mahtkulu maost soolestikku
vdot_maks_veri = m.Param(value=maksu_mass) # l/min, vere mahtkulu maksasse on umbes 1 l / kg maksa. Allikas: C. Eipel, K. Abshagen, and B. Vollmar, “Regulation of hepatic blood flow: The hepatic arterial buffer response revisited,” World Journal of Gastroenterology : WJG, vol. 16, no. 48, p. 6046, Dec. 2010, doi: 10.3748/wjg.v16.i48.6046.
h_magu = m.Param(value=0.0005) # min^-1
h_sool = m.Param(value=0.05) # min^-1

Cet_magu_alg = m_etanool_sisse / mm / V_magu.value # arvutada etanooli algkontsentratsioon kõhus

t_lõpp = 120 # min
npts = 200
m.time = np.linspace(0, t_lõpp, npts) # min

# muutujate defineerimine
Cet_magu = m.Var(Cet_magu_alg) # mol/l, etanooli kontsentratsioon maos
Cet_sool = m.Var(0) # mol/l, etanooli kontsentratsioon soolestikus
Cet_veri = m.Var(0) # mol/l, etanooli kontsentratsioon veres
Cet_maks = m.Var(0) # mol/l, etanooli kontsentratsioon maksas

# mudeli võrrandid
r = m.Intermediate(-Vmax * Cet_maks / (Km + Cet_maks))
ndot_maost = m.Intermediate(h_magu * (Cet_magu - Cet_veri)) # mol/min
ndot_sool = m.Intermediate(h_sool * (Cet_sool - Cet_veri)) # mol/min
m.Equation(Cet_magu.dt() == (-ndot_maost - vdot_magu_sool * Cet_magu) / V_magu) # tegelikult mao maht pidevalt muutub, aga siin lihtsustame ja oletame, et maht ei muutu.
m.Equation(Cet_sool.dt() == (vdot_magu_sool * Cet_magu - ndot_sool) / V_sool)
m.Equation(Cet_veri.dt() == (ndot_maost + ndot_sool - vdot_maks_veri * Cet_veri + vdot_maks_veri * Cet_maks) / V_veri)
m.Equation(Cet_maks.dt() == (vdot_maks_veri * Cet_veri - vdot_maks_veri * Cet_maks + r * V_maks) / V_maks)

m.options.IMODE = 4
m.solve(disp=False)

plt.figure()
plt.plot(m.time, np.asarray(Cet_veri.value) * mm / vere_tihedus)
plt.xlabel('Aeg (min)')
plt.ylabel('Kontsentratsioon (mg/g)')
plt.show()
