'''
Mõned reaktsioonid metüültsüklopentadienüül mangaantrikarbonüüli tootmises.

Viide: R. J. Willey, H. S. Fogler, and M. B. Cutlip, “The integration of process safety into a chemical reaction engineering course: Kinetic modeling of the T2 incident,” Process Safety Progress, vol. 30, no. 1, pp. 39–44, 2011, doi: 10.1002/prs.10431.
'''
import numpy as np
import matplotlib.pyplot as plt

RGAS = 8.314 # ideaalgaasi konstant, J mol^-1 K^-1

def reaktsioon1(t, Cmtpd, Cna):
    '''
    Esimene reaktsioon metüültsüklopentadienüül mangaantrikarbonüüli tootmiseks.

    C6H8 + Na -> C6H7- + Na+ + 0.5 H2

    Parameetrid
    -----------
    t : float
        Temperatuur (degC)
    Cmtpd : float
        Metüültsüklopentadieeni kontsentratsioon (mol dm^-3)
    Cna : float
        Naatriumi kontsentratsioon (mol dm^-3)

    Tagastab
    --------
    r : float
        Reaktsiooni kiirus (mol dm^-3 s^-1)
    Qrxn : float
        Soojusvoog reaktsiooni tõttu (J dm^-3 s^-1)
    '''
    Ea = 40000 # aktivatsioonienergia, J mol^-1 K^-1
    A = 400 / 60 # sagedustegur, dm^3 mol^-1 s^-1
    hrxn = -45400 # reaktsioonientalpia, J mol^-1

    k = A * np.exp(-Ea / RGAS / (t + 273.15)) # t ühik peaks olema siin K, mitte degC
    r = -k * Cmtpd * Cna # mol dm^-3 s^-1
    Qrxn = hrxn * -r # soojusvoog reaktsiooni tõttu, J dm^-3 s^-1
    return r, Qrxn

def reaktsioon2(t, Cdg):
    '''
    Kõrval reaktisoon milles lahusti (diglüüm) laguneb.

    C6H14O3 -> 3 mooli gaasi + jääk

    Parameetrid
    -----------
    t : float
        Temperatuur (degC)
    Cdg : float
        Diglüüm kontsentratsioon (mol dm^-3)

    Tagastab
    --------
    r : float
        Reaktsiooni kiirus (mol dm^-3 s^-1)
    Qrxn : float
        Soojusvoog reaktsiooni tõttu (J dm^-3 s^-1)
    '''
    Ea = 80000 # J mol^-1 K^-1
    A = 1e6 / 60 # dm^3 mol^-1 s^-1
    hrxn = -390000 # reaktsioonientalpia, J mol^-1

    k = A * np.exp(-Ea / RGAS / (t + 273.15)) # t ühik peaks olema siin K, mitte degC
    r = -k * Cdg # mol dm^-3 s^-1
    Qrxn = hrxn * -r # soojusvoog reaktsiooni tõttu, J dm^-3 s^-1
    return r, Qrxn


Cmtpd = 5.1 # mol dm^-3
Cna = 4.3 # mol dm^-3
Cdg = 3 # mol dm^-3
T = np.linspace(0, 350, 100) # degC

r1, Q1 = reaktsioon1(T, Cmtpd, Cna)
r2, Q2 = reaktsioon2(T, Cdg)
Qkokku = Q1 + Q2

plt.figure()
plt.plot(T, -r1, label='metüültsüklopentadieen')
plt.plot(T, -r2, label='diglüüm')
plt.xlabel('Temperatuur (°C)')
plt.ylabel('Reaktsiooni kiirus (mol dm$^{-3}$ s$^{-1}$)')
plt.legend(frameon=False)
plt.show()

plt.figure()
plt.plot(T, -Qkokku)
plt.xlabel('Temperatuur (°C)')
plt.ylabel('Soojuse tekke (J dm^$^{-3}$ s$^{-1}$)')
plt.show()
