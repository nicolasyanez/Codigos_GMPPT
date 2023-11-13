import numpy as np
from matplotlib import pyplot as plt  #NOW LETS MAKE SOME PLOTS.
from pvmismatch import *  #THIS IMPORTS EVERYTHING WE NEED.

plt.style.use("ggplot") #USED FOR HAVE A GRID IN THE PLOT.

pvconst=pvconstants.PVconstants(npts=500) #DEFINITION OF THE NUMBER OF POINTS.
pvcell = pvcell.PVcell(Rs=2.63810001e-04, Rsh= 2.13008168e+01, Isat1_T0=2.74125847e-06, 
                       Isat2_T0=5.51382460e-10, Isc0_T0=10.61 , 
                       aRBD=0.0001036748445065697, bRBD=0.0, VRBD=-5.527260068445654, 
                       nRBD=3.284628553041425, Eg=1.1, alpha_Isc=0.0003551, Tcell=298.15, 
                       Ee=1.0, pvconst=pvconst) #THE PV CELL IS DEFINED WITH ITS CHARACTERISTIC PARAMETERS.
cell_pos = pvmodule.standard_cellpos_pat(20, [2,2,2]) #THE NUMBER OF CELLS AND THE WAY THEY ARE ORDERED ARE DEFINED.
pv_mod = pvmodule.PVmodule(cell_pos=cell_pos, pvcells = pvcell) #THE PHOTOVOLTAIC MODULE IS GENERATED FROM THE DEFINED CELLS.
pv_str = pvstring.PVstring(numberMods=1, pvmods = pv_mod, pvconst = pvconst) #THE PV ARRANGEMENT IS DEFINED AND GENERATED.
pvsys = pvsystem.PVsystem(numberStrs=1, numberMods=1, pvstrs = pv_str) #THE PV SYSTEM IS GENERATED.


def pvsim(params):
    #THIS FUNCTION RECEIVES A PROFILE OF 3 IRRADIANCES, TEMPERATURE AND VOLTAGE.
    #RETURN A DICTIONARY THAT CONTAINS THE POWER, VOLTAGE AND CURRENT CURVES, THE 
    #MAXIMUM POWER AND CURRENT VALUES AND THE CURRENT CORRESPONDING TO THE DELIVERED VOLTAGE.
    G1, G2, G3, T, Vr = params

    #THE CORRECT IRRADIANCE AND TEMPERATURE IS ENTERED FOR EACH PV CELL.
    pvsys.setSuns({0: {0: [(G1/1000., ) * 40, tuple(range(40))]}})
    pvsys.setSuns({0: {0: [(G2/1000., ) * 40, tuple(range(40,80))]}})
    pvsys.setSuns({0: {0: [(G3/1000., ) * 40, tuple(range(80,120))]}})
    pvsys.setTemps(T + 273.15)
    
    P = pvsys.Psys #THE CURRENT CURVE OBTAINED FROM A VOLTAGE SWEEP IS SAVED FROM THE GIVEN IRRADIANCE.
    indP0 = P>=0. #THE INDICES OF THE ELEMENTS THAT CONTAIN VALUES GREATER OR EQUAL TO ZERO ARE SAVED.
    P = P[indP0] #THE POWER CURVE IS SAVED.
    Vvec = pvsys.Vsys[indP0] #THE VOLTAGE CURVE IS SAVED.
    Ivec = pvsys.Isys[indP0] #THE CURRENT CURVE IS SAVED.
    ind = np.where(Vvec<Vr)[0][-1] #THE INDEX IS OBTAINED WHERE THE VOLTAGE IS EQUAL TO THE ONE ENTERED.
    Ir = Ivec[ind] #THE CURRENT CORRESPONDING TO THE ENTERED VOLTAGE IS OBTAINED.
    indm = np.argmax(P) #THE INDEX THAT CONTAINS THE MAXIMUM POWER VALUE IS OBTAINED.
    Pmax = P[indm] #MAXIMUM POWER IS OBTAINED.
    Vmax = Vvec[indm] #MAXIMUM POWER VOLTAGE IS OBTAINED.
    data = {'Ir': Ir, 'Pmax': Pmax, 'Vmax': Vmax, 'Vvec':Vvec, 'Pvec': P, 'Ivec':Ivec} #THE DICTIONARY CONTAINING THE DATA IS CREATED.
    return data

#BUILD ARRAY OF 220 SHADING PATTERNS WITH VALUES BETWEEN 0 AND 1, AND AVOIDING REPEATED CASES.
spc = np.linspace(1,10,10)*0.1
spc_mat = np.zeros((1,3))
for sp1 in spc:
    for sp2 in spc:
        for sp3 in spc:
            spc_mat = np.append(spc_mat, [[sp1, sp2, sp3]], axis = 0)
            
spc_mat = spc_mat[1:]
spc_mat = np.sort(spc_mat, axis = 1)
spc_mat = np.unique(spc_mat, axis = 0)

T = 25

vmax_vec=[] #ARRAY TO SAVE ALL MAXIMUM POWER VOLTAGE.
pmax_vec=[] #ARRAY TO SAVE ALL MAXIMUM POWER.

corriente=[] #AUXILIARY ARRAY TO SAVE THE ACTUAL CURRENT CURVE.
voltaje=np.linspace(0.1,50,128) #ARRAY CONTAINING A VOLTAGE SWEEP BETWEEN 0 AND 50 V.

for k in range(spc_mat.shape[0]):

    spc = spc_mat[k] #ACTUAL IRRADIANCE         
    params = np.append(spc*1000, [T, 1]) #CONCATENATION OF DATA TO ENTER THE pvsim FUNCTION.
    data = pvsim(params.tolist()) #DICTIONATY OBTAINED BY THE pvsim FUNCTION.
    vmax_vec.append(data['Vmax']) #SAVE ACTUAL MAXIMUM POWER VOLTAGE.
    pmax_vec.append(data['Pmax']) #SAVE ACTUAL MAXIMUM POWER.


    for l in voltaje: #VOLTAGE SWEEP TO RE-SAMPLE THE OBTAINED CURVES.
        mx = np.where(data['Vvec']<l)[0][-1] #THE INDEX IS OBTAINED WHERE THE VOLTAGE IS EQUAL TO 'l'.
        I = data['Ivec'][mx] #THE CURRENT CORRESPONDING TO 'l' IS OBTAINED.
        corriente.append(I) #SAVE THE ACTUAL CURRENT INTO ARRAY.
    plt.plot(voltaje,corriente) #PLOT I-V CURVE.
    corriente=[] #CLEAN AUXILIARY ARRAY.

plt.show() #SHOW ALL PLOTS.
