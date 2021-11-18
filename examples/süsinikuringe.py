'''
Praegu üks oluline probleem maailmas on kliima muutmine kasvuhoonegaaside emissioonide tõttu. Peamine kasvuhoonegaas on süsinikdioksiid. CO2 ei jää ainult atmosfääri ja liigub ka ookeanisse ja organismidesse (taimedesse). Lõpuks enamus seda süsinikdioksiidi jõuab tagasi atmosfääri, aga mingi osa jääb ka maa sisse ja moodustab kivimeid või orgaanilisi aineid mullas.

Kui detailidesse süveneda süsinikuringe on päris keeruline, hõlmates palju reaktsioone ja massilevi protsesse. Lisaks, maailm on suur ja tingimused varieeruvad asukohaga. Allolev joonis IPCC 2001. aasta aruandest annab ülevaade süsinikuringest.

Kõike neid kohti saab käsitleda kui reaktoreid selleks, et modelleerida, kuidas CO2 kontsentratsioon atmosfääris muutub ajas. Kuigi süsinikuringe on keeruline, võime lihtsustada mudelit. Võime näiteks käsitleda kõike ookeane kui üks reaktor, kuigi tegelikult ookeanis on mitu kihti, milles tingimused ja protsessid erinevad. Samuti võime mõned osad mudelist välja jätta, mis oluliselt ei mõjuta tulemust, näiteks CO2 vulkaanidest.

Tee mudelit maailma süsinikuringest ja hinda, mis võiks olla CO2 kontsentratsioon atmosfääris 100 aasta pärast siis, kui lõpetatakse kohe fossiilkütuste põletamist.
'''
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

m = GEKKO()

TUNDI_AASTAS = 8760
mm_C = 12.0107 # g/mol
mm_õhk = 0.79 * 2 * 14.0067 + 0.21 * 2 * 15.999 # g/mol
ndot_põletamine = m.Param(value=0 / mm_C) # Pmol C / aasta
k_fotosüntees = m.Param(value=0.11) # mol C taimede kasvule / mol C juba olemas taimedes / aasta. Põhineb andmetel peatükist Allikas: I. C. Prentice et al., “The carbon cycle and atmospheric carbon dioxide,” in Climate change 2001: the scientific basis, International Panel on Climate Change, 2001. Accessed: Nov. 05, 2021. [Online]. Available: https://hal.archives-ouvertes.fr/hal-03333974
k_surm = m.Param(value=0.11) # mol C taimede surmast / mol C olemas taimedes / aasta
ndot_settimine = m.Param(value= 0.5 / mm_C) # Pmol C, mis aastas jäävad mere settesse. Allikas: I. C. Prentice et al., “The carbon cycle and atmospheric carbon dioxide,” in Climate change 2001: the scientific basis, International Panel on Climate Change, 2001. Accessed: Nov. 05, 2021. [Online]. Available: https://hal.archives-ouvertes.fr/hal-03333974
k_lagunemine = m.Param(value=0.04) # mol C tagasi atmosfääri / mol C hetkel mullas / aasta
m_ookean = m.Param(value=1.4e9) # Pg, ookeanite mass. Allikas: D. Avijeet, “Mass of the Oceans,” The Physics Factbook, 1998. https://hypertextbook.com/facts/1998/AvijeetDut.shtml (accessed Nov. 05, 2021).
n_atm = m.Param(value= 5.148e6 / mm_õhk) # Pmol atmosfääris. Allikas: K. E. Trenberth and L. Smith, “The Mass of the Atmosphere: A Constraint on Global Analyses,” Journal of Climate, vol. 18, no. 6, pp. 864–875, Mar. 2005, doi: 10.1175/JCLI-3299.1.
pind_ookean = m.Param(value=360.663e6) # km^2, merede pindala. Allikas: M. J. Costello, A. Cheung, and N. De Hauwere, “Surface Area and the Seabed Area, Volume, Depth, Slope, and Topographic Variation for the World’s Seas, Oceans, and Countries,” Environ. Sci. Technol., vol. 44, no. 23, pp. 8821–8828, Dec. 2010, doi: 10.1021/es1012752.
k = m.Param(value=0.0001) # km h^-1, Maa keskmine gaasi läbikandekiirus atmosfääri ja ookeani vahel. Allikas: F. Thomas, C. Perigaud, L. Merlivat, and J.-F. Minster, “World-scale monthly mapping of the C02 ocean-atmosphere gas-transfer coefficient,” Philosophical Transactions of the Royal Society of London. Series A, Mathematical and Physical Sciences, vol. 325, no. 1583, pp. 71–83, 1988.
sol = m.Param(value=32.93e-6) # Pmol km^-3 bar^-1, CO2 lahustuvus merevees 20 kraadi juures. Allikas: F. Thomas, C. Perigaud, L. Merlivat, and J.-F. Minster, “World-scale monthly mapping of the C02 ocean-atmosphere gas-transfer coefficient,” Philosophical Transactions of the Royal Society of London. Series A, Mathematical and Physical Sciences, vol. 325, no. 1583, pp. 71–83, 1988.

t_lõpp = 100 # aastat
npts = 200
m.time = np.linspace(0, t_lõpp, npts)

CO2_2020 = 413 # ppm
CO2_ookean_nüüd = 0.0024 # mol/kg, Allikas: R. Zeebe and D. Wolf-Gladrow, “Carbon dioxide, dissolved (Ocean),” in Encyclopedia of paleoclimatology and ancient environments / ed. by Vivien Gornitz, Dordrecht: Springer, 2009, pp. 123–127. Accessed: Nov. 05, 2021. [Online]. Available: https://epic.awi.de/id/eprint/11770/
CCO2_atm = m.Var(CO2_2020) # ppm
CCO2_ookean = m.Var(CO2_ookean_nüüd) # mol/kg
nC_taimed = m.Var(500/mm_C) # Pmol, taimede süsiniku kogus
nC_muld = m.Var(1350/mm_C) # Pmol (10^15 mol), Allikas: I. C. Prentice et al., “The carbon cycle and atmospheric carbon dioxide,” in Climate change 2001: the scientific basis, International Panel on Climate Change, 2001. Accessed: Nov. 05, 2021. [Online]. Available: https://hal.archives-ouvertes.fr/hal-03333974

PCO2_ookean = m.Intermediate(165290 * CCO2_ookean) # CO2 osarõhk ookeani pinnas
ndot_lahustumine = m.Intermediate(k * TUNDI_AASTAS * sol * pind_ookean * (CCO2_atm - PCO2_ookean) / 1e6)
m.Equation(CCO2_atm.dt() == (ndot_põletamine - k_fotosüntees * nC_taimed + k_lagunemine * nC_muld - ndot_lahustumine) / n_atm * 1e6)
m.Equation(CCO2_ookean.dt() == (ndot_lahustumine - ndot_settimine) / m_ookean * 1000)
m.Equation(nC_taimed.dt() == k_fotosüntees * nC_taimed - k_surm * nC_taimed)
m.Equation(nC_muld.dt() == k_surm * nC_taimed - k_lagunemine * nC_muld)

m.options.IMODE = 4
m.solve(disp=False)

plt.figure()
plt.plot(m.time, CCO2_atm.value)
plt.xlabel('Aeg (aasta)')
plt.ylabel('CO$_2$ atmosfääris (ppm)')
plt.show()
