import math

def ph(h_plus):
    return -math.log10(h_plus)

def residence_time(Vriv,Criv,Vocean,Cocean):
    return (Vocean*Cocean)/(Vriv*Criv)

def h_plus(ph):
    return math.pow(10,-ph)

def HCO3(ALK,DIC):
    return 2*DIC-ALK

def CO3(ALK,DIC):
    return ALK-DIC

def pCO2(k0,k1,k2,ALK,DIC):
    return (k2*math.pow(2*(DIC-ALK),2))/(k0*k1*(ALK-DIC))

def saturation_state(Ca,CO3,ksp):
    return (Ca*CO3)/ksp

def ph2(k2,ALK,DIC):
    return k2*(2*DIC-ALK)/(ALK-DIC)

def bio_removal(Vriv,Criv,Vmix,Cdeep,Csurf):
    return (Vriv*Criv+Vmix*Cdeep-Vmix*Csurf)/(Vriv*Criv+Vmix*Cdeep)

def burial(Vriv,Criv,Vmix,Cdeep,Csurf):
    return (Vriv*Criv)/(Vriv*Criv+Vmix*Cdeep-Vmix*Csurf)

def light_tensity(z,a,I0):
    return I0*math.exp(-a*z)

def Michaelis_Menten(Vmax,Cs,Km):
    return (Vmax*Cs)/(Km+Cs)

def Biological_pump(Cdeep,Csurf):
    return (Cdeep-Csurf)/Cdeep

def stokes_sink(g,viscosity,radius,density_diff):
    return 2*g*math.pow(radius,2)*(density_diff)/(9*viscosity)

def Martin(initial,z):
    return initial*math.pow(z/100,-0.858)

def GAS_LAW(p,v,n,t,r):
    return p*v == n*t*r

def Henry(Kh,c):
    return Kh*c

def AOU(T,S,O2):
    return math.exp(177.7888 + 255.5907 /(T/100) + 146.4813* math.log(T/100,math.e)-22.204 * (T/ S*( 0.037362 +
            (T/100)* (0.016504-0.0020564 * T/100)))* 1000/22.392)-O2

def multiplicative_nutrient(Vmax,Fe,Si,Kmfe,Kmsi):
    return Vmax*(Fe/(Kmfe+Fe))*(Si/(Kmsi+Si))

def Vmix(C14surf,C14deep,Lambda,Vdeep):
    return Lambda*Vdeep/(C14surf/C14deep-1)

def Activity(Lambda,n):
    return Lambda*n

def radioactive_decay(N0,Lambda,t):
    return N0*math.exp(-Lambda*t)

def radioactive_growth(Np,Lambda,t):
    return Np*(math.exp(Lambda*t)-1)

def mean_life(Lambda):
    return 1/Lambda

def half_life(Lambda):
    return math.log(2,math.e)/Lambda

def free_metal(Ctotal,logbeta_list,C_list,EDTA,Rinorganic):
    count = 1
    count += C_list[0]*logbeta_list[0]
    count += math.pow(C_list[0],2)*logbeta_list[1]
    count += C_list[1]*math.pow(10,logbeta_list[2])
    count += C_list[2]*math.pow(10,logbeta_list[3])
    Y4 = EDTA/count
    result = Ctotal/(math.pow(10,Rinorganic)+math.pow(10,logbeta_list[4])*Y4)
    ratio = result/Ctotal*100
    return result,ratio