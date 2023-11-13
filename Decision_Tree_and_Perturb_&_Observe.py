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

#************* Desition Tree **************/

def DT(vPV,iPV):#RECIEVE VOLTAGE AND CURRENT IN PV SYSTEM AND RETURN A VOLTAGE REFERENCE.

    clase=np.array([40,20,10])
    leaf=0
    state=0
    P=vPV*iPV
    if(vPV >= 30.0):
        state=0
    if(vPV >= 15.0 and vPV < 30):
        state=1
    if(vPV < 15.0):
        state=2
    match state:
    #//////////////////////////////////////////////////////////////////////////
    #///////V40
        case 0:
            if (P <= 132.77407455444336): 
                if (vPV <= 39.13030433654785):
                    if (P <= 52.92519760131836):
                        if (P <= 35.72148895263672):
                            leaf=2
                        else:
                            if (P <= 39.71645927429199):
                                if (vPV <= 37.23977851867676):
                                    leaf=1
                                else:
                                    if (vPV <= 37.81496047973633):
                                        leaf=0
                                    else:
                                        if (vPV <= 38.04804611206055):
                                            leaf=1
                                        else:
                                            if (vPV <= 38.3769474029541):
                                                if (vPV <= 38.20389938354492):
                                                    if (vPV <= 38.126041412353516):
                                                        leaf=0
                                                    else:
                                                        leaf=1
                                                else:
                                                    leaf=0
                                            else:
                                                if (P <= 39.50600242614746):
                                                    leaf=1
                                                else:
                                                    leaf=0
                            else:
                                leaf=1
                    else:
                        if (vPV <= 38.554697036743164):
                            if (vPV <= 38.01014518737793):
                                leaf=2
                            else:
                                if (vPV <= 38.190975189208984):
                                    leaf=0
                                else:
                                    leaf=2
                        else:
                            if (P <= 99.13892364501953):
                                if (vPV <= 38.65224075317383):
                                    leaf=1
                                else:
                                    leaf=0
                            else:
                                leaf=2         
                else:
                    if (P <= 80.74540328979492):
                        if (P <= 80.72513961791992):
                            if (vPV <= 39.32780838012695):
                                if (vPV <= 39.2650203704834):
                                    leaf=1
                                else:
                                    leaf=0
                            else:
                                leaf=1
                        else:
                            leaf=0
                    else:
                        leaf=1
            else:
                if (vPV <= 40.26659393310547):
                    if (vPV <= 40.06301307678223):
                        leaf=2
                    else:
                        if (P <= 185.6413116455078):
                            leaf=1
                        else:
                            leaf=2
                else:
                    if (P <= 208.67066192626953): 
                        leaf=1
                    else:
                        leaf=2
    #//////////////////////////////////////////////////////////////////////////
    #///////V10
        case 2:
            if(vPV > 12.0):
                leaf=1
            else:
                leaf=2
    return clase[leaf]


def po(vact,vpas,ipas,step):
    #FUNCTION THAT PLAYS THE P&O
    params = np.append(current_irr,[T,vact])
    data = pvsim(params.tolist()) #PV DATA IS OBTAINED
    I_pv=data['Ir'] #PV CURRENT AT VACT INSTANT IS EXTRACTED
    pact=vact*I_pv #ACTUAL PV POWER IS CALCULATED
    ppas=vpas*ipas #PASS PV POWER IS CALCULATED
    #P&O
    if pact > ppas:                               #IF P_ACT IS GREATER THAN P_PAS
        if vact > vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #STEP UP REFERENCE
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #STEP DOWN REFERENCE
            return vact,vpas,ipas
    else:                                            #IF P_ACT IS NOT GREATER THAN P_PAS
        if vact < vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #STEP UP REFERENCE
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #STEP DOWN REFERENCE
            return vact,vpas,ipas


current_irr=irr[0,:] #CURRENT IRRADIANCE

valoresV=np.linspace(5,40,8)
auxP=np.zeros((8))

step=1/156
maq=0
gmpptaux=0
poaux=0

vact=40 #INITIALIZE P&O PARAMETERS
vpas=0
ipas=0
step=0.2
vact,vpas,ipas=po(vact,vpas,ipas,step) #INITIALIZE P&O

m=np.ones(8)*10 #ARRAY TO WAIT FOR STABILIZATION
stab_m=0 #STABILIZATION COUNTER


for instance in range(156):

    if(instance<52): #FUNCTION TO CHANGE CURRENT IRRADIANCE
        current_irr=irr[0,:]
    if(instance>=52 and instance<104):
        current_irr=irr[1,:]
    if(instance>=104):
        current_irr=irr[2,:]
    if(instance==52 or instance==104): #AT A FIXED TIME THE STATE MACHINE CHANGES
        maq=1

    match maq:
        case 0: #P&O STAGE
            vact,vpas,ipas=po(vact,vpas,ipas,step)
            vvec[iter]=vact
            ivec[iter2]=ipas
            iter=iter+1
            iter2=iter2+1
        case 1: #DECISION TREE STAGE
            vact,vpas,ipas=po(vact,vpas,ipas,step) #START WAITING TO P&O STABILIZATION
            vvec[iter]=vact
            ivec[iter2]=ipas
            iter=iter+1
            iter2=iter2+1
            m[stab_m]=(vact-vpas)/step
            stab_m+=1
            if(stab_m==8):
                stab_m=0 
            suma_m=0
            suma_m=np.sum(m)
            if(suma_m == 0): #IT MEANS THAT IF P&O IS STABILIZED
                m=np.ones(8)*10
                vact=DT(vact,ipas) #RUN DECISION TREE
                stab_m=0
                maq=0

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


# GRAPHIC 1
axes[0].step(x, vvec*ivec, where='post')
axes[0].step(x, pmax, where='post', color='red')
axes[0].set_title('Potencia P(t)')
axes[0].grid(axis='x', color='0.95')
axes[0].grid(axis='y', color='0.95')
axes[0].set_xlim([0, 6])     # X AXIS LIMITS
axes[0].set_ylim([0, 400])   # Y AXIS LIMITS


# GRAPHIC 2
axes[1].step(x, vvec, where='post')
axes[1].set_title('Voltaje V(t)')
axes[1].grid(axis='x', color='0.95')
axes[1].grid(axis='y', color='0.95')
axes[1].set_xlim([0, 6])    # X AXIS LIMITS
axes[1].set_ylim([0, 50])   # Y AXIS LIMITS


# GRAPHIC 3
axes[2].step(x, ivec, where='post')
axes[2].set_title('Corriente I(t)')
axes[2].grid(axis='x', color='0.95')
axes[2].grid(axis='y', color='0.95')
axes[2].set_xlim([0, 6])    # X AXIS LIMITS
axes[2].set_ylim([0, 15])   # Y AXIS LIMITS

#plt.savefig("dtpo.pdf")

plt.show()
