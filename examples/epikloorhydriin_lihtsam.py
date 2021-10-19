'''
Epikloorhüdriin on oluline lähteaine epoksiidvaikude, elastomeeride ja sünteetilise glütseriini tootmiseks. Seda praegu toodetakse enamasti propüleenist ja kloorgaasist, aga on võimalik ka toota glütseroolist. Tootmine glütseroolist on võibolla eelistatud kuna ei ole vaja kasutada kloorgaasi. Esimeses sammus toimub reaktsioon kloorhappe ja glütserooli vahel katalüsaatori juuresolekul selleks, et toota 2,3-dikloro-1-propanooli (DKP). Seejärel viiakse läbi reaktsiooni DKP ja naatrium hüdroksiidi vahel epikloorhüdriini sünteesimiseks.

C3H6OCl2 + NaOH -> C3H5OCl + H2O + NaCl

Tahetakse läbi viia selle reaktsiooni 1-liitrilises perioodilises reaktoris. Kui kaua aega segu peaks olema reaktoris selleks, et saavutada 90% konversiooniastet?

DKP algkontsentratsioon on 0,08 mol/l ja naatrium hüdroksiidi algkontsentratsioon on 0,08 mol/l. Reaktorit hoitakse 60 kraadi juures. DKP suhtes on teist järku reaktsioon ja NaOH suhtes esimest järku. Reaktsiooni kineetilised parameetrid:
A = 1,61*10^25 l^3 mol^-2 s^-1
Ea = 150 kJ/mol

Viide: J. S. Zhang, Y. C. Lu, Q. R. Jin, K. Wang, and G. S. Luo, “Determination of kinetic parameters of dehydrochlorination of dichloropropanol in a microreactor,” Chemical engineering journal, vol. 203, pp. 142–147, 2012.
'''
import numpy as np

RGAS = 8.314462618 # J mol^-1 K^-1

def kiiruskonstant(T):
    A = 1.61e25 # l^3 mol^-2 s^-1
    Ea = 150 * 1000 # J mol^-1
    k = A * np.exp(-Ea / RGAS / (T + 273.15))
    return k # l^3 mol^-2 s^-1

C_dkp_alg = 0.08 # mol l^-1
C_oh_alg = 0.08 # mol l^-1
V = 1 # l
T = 60 # Celcius
X_dkp = 0.9 # konversiooniaste

k = kiiruskonstant(T)
print('kiiruskonstant=', k)

'''
Alustame perioodilise reaktori põhivõrrandiga ja lisame reaktsiooni kiirusevalemi:
    t = -n_dkp_alg \\integraal dX_dkp / (V * (-k * C_dkp_alg^2 * (1 - X_dkp)^2 * (C_oh_alg - C_dkp_alg * X_dkp)))

Kuna aine A moolide arv alguses on lihtsalt maht korrutatud kontsentratsiooniga,
saame kirjutada valemit järgmiselt:
    t = -V * C_dkp_alg \\integraal dX_dkp / (V * (-k * C_dkp_alg^2 * (1 - X_dkp)^2 * (C_oh_alg - C_dkp_alg * X_dkp)))

Võime kõik konstandid integraalist välja tuua ja osa neist taanduvad välja:
    t = 1 / (k * C_dkp_alg) \\integraal dX_dkp / ((1 - X_dkp)^2 * (C_oh_alg - C_dkp_alg * X_dkp))

Kui mõlemad algkontsentratsioonid on samad, siis võime valemit lihtsustada.
    t = 1 / (k * C_dkp_alg^2) \\integraal dX_dkp / (1 - X_dkp)^3

Siis saame integreerida, et saada allolevat valemit.
'''
t = 1 / k / C_dkp_alg**2 * (1 / 2 / (1 - X_dkp)**2 - 1 / 2)
print('Vajalik aeg= {:.2f} s'.format(t))
