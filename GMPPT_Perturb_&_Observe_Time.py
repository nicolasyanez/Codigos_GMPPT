import numpy as np
from matplotlib import pyplot as plt  # now lets make some plots
from pvmismatch import *  # this imports everything we need

plt.style.use("ggplot")


pvconst=pvconstants.PVconstants(npts=2000)
pvcell = pvcell.PVcell(Rs=2.63810001e-04, Rsh= 2.13008168e+01, Isat1_T0=2.74125847e-06, 
                       Isat2_T0=5.51382460e-10, Isc0_T0=10.61 , 
                       aRBD=0.0001036748445065697, bRBD=0.0, VRBD=-5.527260068445654, 
                       nRBD=3.284628553041425, Eg=1.1, alpha_Isc=0.0003551, Tcell=298.15, 
                       Ee=1.0, pvconst=pvconst)
cell_pos = pvmodule.standard_cellpos_pat(20, [2,2,2])
pv_mod = pvmodule.PVmodule(cell_pos=cell_pos, pvcells = pvcell)
pv_str = pvstring.PVstring(numberMods=1, pvmods = pv_mod, pvconst = pvconst)
pvsys = pvsystem.PVsystem(numberStrs=1, numberMods=1, pvstrs = pv_str)


def pvsim(params):
    G1, G2, G3, T, Vr = params
    
    pvsys.setSuns({0: {0: [(G1/1000., ) * 40, tuple(range(40))]}})
    pvsys.setSuns({0: {0: [(G2/1000., ) * 40, tuple(range(40,80))]}})
    pvsys.setSuns({0: {0: [(G3/1000., ) * 40, tuple(range(80,120))]}})
    pvsys.setTemps(T + 273.15)
    P = pvsys.Psys
    indP0 = P>=0.
    P = P[indP0]
    Vvec = pvsys.Vsys[indP0]
    Ivec = pvsys.Isys[indP0]
    ind = np.where(Vvec<Vr)[0][-1]
    Ir = Ivec[ind]
    indm = np.argmax(P)
    Pmax = P[indm]
    Vmax = Vvec[indm]
    data = {'Ir': Ir, 'Pmax': Pmax, 'Vmax': Vmax, 'Vvec':Vvec, 'Pvec': P, 'Ivec':Ivec}
    return data

#ARRAYS ARE CREATED TO SAVE THE POWER, VOLTAGE AND CURRENT CURVES.
vvec=np.ones(156)
ivec=np.ones(156)
pmax=np.ones(156)
pvec=np.ones(156)
iter=0
iter2=0

v=44.2
irr =np.array([[1000,1000,1000],[1000,500,200],[100,100,100]]) #IRRADIANCE PROFILES USED.
T=25.0


def po(vact,vpas,ipas,step):
    #FUNCTION THAT PLAYS THE P&O
    
    params = np.append(current_irr,[T,vact])
    data = pvsim(params.tolist())
    I_pv=data['Ir']
    pact=vact*I_pv
    ppas=vpas*ipas

    #P&O
    if pact > ppas:                               #IF P_ACT IS GREATER THAN P_PAS
        if vact > vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #UP
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #DOWN
            return vact,vpas,ipas
    else:                                            #IF P_ACT IS NOT GREATER THAN P_PAS
        if vact < vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #UP
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #DOWN
            return vact,vpas,ipas


current_irr=irr[0,:]#CURRENT IRRADIANCE

valoresV=np.linspace(5,40,8)
auxP=np.zeros((8))

step=1/156
maq=0
gmpptaux=0
poaux=0

vact=40
vpas=0
ipas=0
step=0.2
vact,vpas,ipas=po(vact,vpas,ipas,step)#INITIALIZE P&O


for instance in range(156):

    if(instance<52):
        current_irr=irr[0,:]
    if(instance>=52 and instance<104):
        current_irr=irr[1,:]
    if(instance>=104):
        current_irr=irr[2,:] #CURRENT IRRADIANCE PROFILE

    if(instance == 78): #WHEN A COUNT IS COMPLETED, IT GOES TO THE SWEEP STAGE
        maq=1

    match maq:
        case 0: #P&O STAGE
            vact,vpas,ipas=po(vact,vpas,ipas,step)
            vvec[iter]=vact
            ivec[iter2]=ipas
            iter=iter+1
            iter2=iter2+1
        case 1:#VOLTAGE SWEEP STAGE
            vvec[iter]=valoresV[gmpptaux]
            parametros = np.append(current_irr,[T,valoresV[gmpptaux]])
            datos = pvsim(parametros.tolist())
            ivec[iter2]=datos['Ir']
            auxP[gmpptaux]=vvec[iter]*ivec[iter2]
            iter=iter+1
            iter2=iter2+1
            gmpptaux=gmpptaux+1
            if(gmpptaux > 7):
                maq=0
                gmpptaux=0
                vact=valoresV[np.argmax(auxP)]

#FUNTION TO PLOT MAXIMUM POWER IN TIME.
for z in range(156):
    if(z<52):
        pmax[z]=380
    if(z>=52 and z<104):
        pmax[z]=135
    if(z>=104):
        pmax[z]=33

x=np.arange(0,156*(1/26),1/26)

fig, axes = plt.subplots(nrows=3, ncols=1)#, figsize=(8, 40))


# Gráfico 1
axes[0].step(x, vvec*ivec, where='post')
axes[0].step(x, pmax, where='post', color='red')
axes[0].set_title('Potencia P(t)')
axes[0].grid(axis='x', color='0.95')
axes[0].grid(axis='y', color='0.95')
axes[0].set_xlim([0, 6])     # Límites del eje x
axes[0].set_ylim([0, 400])   # Límites del eje y


# Gráfico 2
axes[1].step(x, vvec, where='post')
axes[1].set_title('Voltaje V(t)')
axes[1].grid(axis='x', color='0.95')
axes[1].grid(axis='y', color='0.95')
axes[1].set_xlim([0, 6])    # Límites del eje x
axes[1].set_ylim([0, 50])   # Límites del eje y


# Gráfico 3
axes[2].step(x, ivec, where='post')
axes[2].set_title('Corriente I(t)')
axes[2].grid(axis='x', color='0.95')
axes[2].grid(axis='y', color='0.95')
axes[2].set_xlim([0, 6])    # Límites del eje x
axes[2].set_ylim([0, 15])   # Límites del eje y


#plt.savefig("gmpptpo2.pdf")
plt.show()
