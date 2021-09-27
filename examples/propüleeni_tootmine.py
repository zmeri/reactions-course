'''
Praegu puhastatakse reovett enamasti tehnoloogiatega, mis nõuavad energiat ja kemikaale.
Kuid on pakutud, et oleks võimalik toota väärtuslikud kemikaale reoveest. Näiteks, on
leitud baktereid, mis toodavad looduslikku polümeeri polühüdroksübutüraati. Kuumutades
polümeeri on võimalik toota muid kemikaale lagunemise teel. Lagunemises peamine produkt
on 3-hüdroksübutaanhape ja see omakorda võib edasi reageerida moodustamiseks propüleeni,
mis on väärtuslik vaheprodukt keemiatööstuses.

C4H8O3 -> C3H6 + H2O + CO2

Kui palju propüleeni oleks võimalik toota, kui töödelda 1000 kuupmeetrit lahust päevas
100 liitrises pidev segureaktoris? 3-hüdroksübutaanhappe kontsentratsioon alglahuses on 2,2 mol/l
ja reaktsioon toimub 250 kraadi juures.

Kineetilised parameetrid:
ln(A) = 29,2
Ea = 126,8 kJ/mol

Viide: Y. Li and T. J. Strathmann, “Kinetics and mechanism for hydrothermal conversion
of polyhydroxybutyrate (PHB) for wastewater valorization,” Green Chemistry, vol. 21,
no. 20, pp. 5586–5597, 2019, doi: 10.1039/C9GC02507C.
'''
import numpy as np

RGAS = 8.314462618 # J mol^-1 K^-1

def kiiruskonstant(T):
    A = np.exp(29.2)
    Ea = 126.8 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15)) # s^-1
    return k

Chbh_alg = 2.2 # mol l^-1
V = 100 # l
T = 250 # Celcius
vdot = 1000 # m^-3 päev^-1
vdot = vdot * 1000 / 24 / 3600 # l s^-1

k = kiiruskonstant(T)

'''
Kasutades segureaktori põhivõrrandit (mis tuleb moolbilansist) koostame järgmist valemit:
    V = ndot_hbh_alg * Xhbh / (k * Chbh_alg * (1 - Xhbh))

Kuna 3-hüdroksübutaanhappe moolkulu alguses on lihtsalt mahtkulu korrutatud kontsentratsiooniga,
saame kirjutada valemit järgmiselt:
    V = vdot * Chbh_alg * Xhbh / (k * Chbh_alg * (1 - Xhbh))

Lahendades Xhbh jaoks saame allolevat valemit.
'''
Xhbh = V * k * Chbh_alg / (vdot * Chbh_alg + V * k * Chbh_alg)
print('3-hüdroksübutaanhappe konversiooniaste= {:.2f}%'.format(Xhbh * 100))

mm_prop = 42.08 # propüleeni molaar mass, g mol^-1
mdot_prop = vdot * Chbh_alg * Xhbh * mm_prop / 1000 # kg s^-1
print('Toodetud propüleen= {} kg/s'.format(mdot_prop))
print('Toodetud propüleen= {} kg/päev'.format(mdot_prop * 3600 * 24))
