#Kint calculation values as determined in Nguyen D, Mayne L, Phillips MC, Walter Englander S. Reference Parameters for
# Protein Hydrogen Exchange Rates. J Am Soc Mass Spectrom. 2018;29(9):1936-1939. doi:10.1007/s13361-018-2021-z
from math import exp, log10

#Log10kex(X) values
params = {
    'A' : [0.00, 0.00, 0.00, 0.00],
    'R' : [-0.59, -0.32, 0.08, 0.22],
    'N' : [-0.58, -0.13, 0.49, 0.32],
    'C' : [-0.54, -0.46, 0.62, 0.55],
    'G' : [-0.22, 0.22, -0.03, 0.17],
    'Q' : [-0.47, -0.27, 0.06, 0.20],
    'I' : [-0.91, -0.59, -0.73, -0.23],
    'L' : [-0.57, -0.13, -0.58, -0.21],
    'K' : [-0.56, -0.29, -0.04, 0.12],
    'M' : [-0.64, -0.28, -0.01, 0.11],
    'F' : [-0.52, -0.43, -0.24, 0.06],
    'P' : ['NA', -0.19, 'NA', -0.24],
    'S' : [-0.44, -0.39, 0.37, 0.30],
    'T' : [-0.79, -0.47, -0.07, 0.20],
    'W' : [-0.40, -0.44, -0.41, -0.11],
    'Y' : [-0.41, -0.37, -0.27, 0.05],
    'V' : [-0.74, -0.30, -0.70, -0.14]
}
Nterm_acid = -1.32
Nterm_base = 1.62

Cterm_acid=0.96
Cterm_base=-1.80

#Activation energies, used to calculate HX (Kx) rates at different temperatures
Ea = 14.000 #Ea(Ka) = 14 kcal/mol
Eb = 17.000 #Eb(Kb) = 17 kcal/mol
Ew = 19.000 #Ew(Kw) = 19 kcal/mol

#Residue specific activation energies
Ea_Asp = 1.000 #kcal/mol
Ea_Glu = 1.083 #kcal/mol
Ea_His = 7.500 #kcal/mol

#Reference Ala-Ala exchange rates at 20C
refKa = 10**1.62 #Ka, M-1min-1
refKb = 10**10.18 #Kb, M-1min-1
refKw = 10**-1.5 #Kw, M-1min-1

#Temperature dependence of HX rates
R = 0.001987 #kcal/molK, universal gas constant
refKwD20 = 15.05 #molar ionisation constant for D20 at 20C

#Calc Kint
def calc_kint(res1, res2, temperature, pD):

    Ka_tmod = refKa * exp(-Ea * ((1 / temperature) - (1 / 293)) / R)
    Kb_tmod = refKb * exp(-Eb * ((1 / temperature) - (1 / 293)) / R)
    Kw_tmod = refKw * exp(-Ew * ((1 / temperature) - (1 / 293)) / R)

    #logk(acid) = logkA,ref + logAL + logAR - pD.
    logkAcid = Ka_tmod * 10**(params[res1][0] + params[res2][1]) * 10**-pD
    logkBase = Kb_tmod * 10**(params[res1][2] + params[res2][3]) * 10**(pD-refKwD20)
    logkWater = Kw_tmod * 10**(params[res1][2] + params[res2][3])

    return(logkAcid + logkBase + logkWater)

def calc_kint_Alt(res1, res2, temperature, pD, pos, length):

    His_pKc = -log10((10**-7.42)*exp(-1*Ea_His*(((1/temperature)-(1/278))/R)))
    Glu_pKc = -log10((10**-4.93)*exp(-1*Ea_Glu*(((1/temperature)-(1/278))/R)))
    Asp_pKc = -log10((10**-4.48)*exp(-1*Ea_Asp*(((1/temperature)-(1/278))/R)))

    if res2 == 'H':
        Acid_Kex_R = log10(((10**(-0.51-pD))/((10**-His_pKc)+(10**-pD))) + \
                           ((10**(-His_pKc))/((10**-His_pKc)+(10**-pD)))) #cannot acid catalyse with His, only His+?
    elif res2 == 'E':
        Acid_Kex_R = log10((10 ** (-0.27 - pD)) / ((10 ** -Glu_pKc) + (10 ** -pD)) + \
                    (10 ** (0.31-Glu_pKc)) / ((10 ** -Glu_pKc) + (10 ** -pD)))
    elif res2 == 'D':
        Acid_Kex_R = log10((10 ** (-0.12 - pD)) / ((10 ** -Asp_pKc) + (10 ** -pD)) + \
                    (10 ** (0.58-Asp_pKc)) / ((10 ** -Asp_pKc) + (10 ** -pD)))
    else:
        Acid_Kex_R = params[res2][1]

    if res1 == 'H':
        Acid_Kex_L = log10((10**(-0.8-pD))/((10**-His_pKc)+(10**-pD)) + \
                    (10**(-His_pKc))/((10**-His_pKc)+(10**-pD)))
    elif res1 == 'E':
        Acid_Kex_L = log10((10 ** (-0.6 - pD)) / ((10 ** -Glu_pKc) + (10 ** -pD)) + \
                    (10 ** (-0.9-Glu_pKc)) / ((10 ** -Glu_pKc) + (10 ** -pD)))
    elif res1 == 'D':
        Acid_Kex_L = log10((10 ** (-0.9 - pD)) / ((10 ** -Asp_pKc) + (10 ** -pD)) + \
                    (10 ** (0.9-Asp_pKc)) / ((10 ** -Asp_pKc) + (10 ** -pD)))
    else:
        Acid_Kex_L = params[res1][0]

    if res2 == 'H':
        Base_Kex_R = log10((10 ** (0.83 - pD)) / ((10 ** -His_pKc) + (10 ** -pD)) + \
                           (10 ** (0.14-His_pKc)) / ((10 ** -His_pKc) + (10 ** -pD)))  # cannot acid catalyse with His, only His+?
    elif res2 == 'E':
        Base_Kex_R = log10((10 ** (-0.39 - pD)) / ((10 ** -Glu_pKc) + (10 ** -pD)) + \
                           (10 ** (-0.15 - Glu_pKc)) / ((10 ** -Glu_pKc) + (10 ** -pD)))
    elif res2 == 'D':
        Base_Kex_R = log10((10 ** (0.6 - pD)) / ((10 ** -Asp_pKc) + (10 ** -pD)) + \
                           (10 ** (-0.18 - Asp_pKc)) / ((10 ** -Asp_pKc) + (10 ** -pD)))
    else:
        Base_Kex_R = params[res2][3]

    if res1 == 'H':
        Base_Kex_L = log10((10 ** (0.8 - pD)) / ((10 ** -His_pKc) + (10 ** -pD)) + \
                           (10 ** (-0.1-His_pKc)) / ((10 ** -His_pKc) + (10 ** -pD)))
    elif res1 == 'E':
        Base_Kex_L = log10((10 ** (0.24 - pD)) / ((10 ** -Glu_pKc) + (10 ** -pD)) + \
                           (10 ** (-0.11 - Glu_pKc)) / ((10 ** -Glu_pKc) + (10 ** -pD)))
    elif res1 == 'D':
        Base_Kex_L = log10((10 ** (0.69 - pD)) / ((10 ** -Asp_pKc) + (10 ** -pD)) + \
                           (10 ** (0.1 - Asp_pKc)) / ((10 ** -Asp_pKc) + (10 ** -pD)))
    else:
        Base_Kex_L = params[res1][2]

    if pos == 1:
        Acid_Kex_R = Acid_Kex_R + Nterm_acid
        Base_Kex_R = Base_Kex_R + Nterm_base
    elif pos == length-1:
        Acid_Kex_L = Acid_Kex_L + Cterm_acid
        Base_Kex_L = Base_Kex_L + Cterm_base

    Ka_tmod = refKa * exp(-Ea * ((1 / temperature) - (1 / 293)) / R)
    Kb_tmod = refKb * exp(-Eb * ((1 / temperature) - (1 / 293)) / R)
    Kw_tmod = refKw * exp(-Ew * ((1 / temperature) - (1 / 293)) / R)

    if pos == 0:
        return -1
    elif res1 == 'P':
        return -1
    else:
        logkAcid = Ka_tmod * 10 ** (Acid_Kex_L + Acid_Kex_R) * 10 ** -pD
        logkBase = Kb_tmod * 10 ** (Base_Kex_L + Base_Kex_R) * 10 ** (pD - refKwD20)
        logkWater = Kw_tmod * 10 ** (Base_Kex_L + Base_Kex_R)
        return(logkAcid + logkBase + logkWater)

