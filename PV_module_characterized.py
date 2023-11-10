import numpy as np
from matplotlib import pyplot as plt  # now lets make some plots
from pvmismatch import *  # this imports everything we need

plt.style.use("ggplot")

pvconst=pvconstants.PVconstants(npts=500) #DEFINICIÓN DE LA CANTIDAD DE PUNTOS
pvcell = pvcell.PVcell(Rs=2.63810001e-04, Rsh= 2.13008168e+01, Isat1_T0=2.74125847e-06, 
                       Isat2_T0=5.51382460e-10, Isc0_T0=10.61 , 
                       aRBD=0.0001036748445065697, bRBD=0.0, VRBD=-5.527260068445654, 
                       nRBD=3.284628553041425, Eg=1.1, alpha_Isc=0.0003551, Tcell=298.15, 
                       Ee=1.0, pvconst=pvconst) #SE DEFINE LA CELDA FV CON SUS PARÁMETROS CARACTERÍSTICOS
cell_pos = pvmodule.standard_cellpos_pat(20, [2,2,2]) #SE DEFINE LA CANTIDAD DE CELDAS Y DE SU ORDENAMIENTO
pv_mod = pvmodule.PVmodule(cell_pos=cell_pos, pvcells = pvcell) #SE GENERA EL MÓDULO FOTOVOLTAICO A PARTIR DE LAS CELDAS DEFINIDAS
pv_str = pvstring.PVstring(numberMods=1, pvmods = pv_mod, pvconst = pvconst) #SE DEFINE Y GENERA EL ARREGLO FV
pvsys = pvsystem.PVsystem(numberStrs=1, numberMods=1, pvstrs = pv_str) #SE GENERA EL SISTEMA FV


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

#Build array of shading patterns
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

vmax_vec=[]
pmax_vec=[]

corriente=[]
voltaje=np.linspace(0.1,50,128)

for k in range(spc_mat.shape[0]):

    spc = spc_mat[k]          
    params = np.append(spc*1000, [T, 1])
    data = pvsim(params.tolist())
    vmax_vec.append(data['Vmax'])
    pmax_vec.append(data['Pmax'])


    for l in voltaje:
        mx = np.where(data['Vvec']<l)[0][-1]
        I = data['Ivec'][mx]
        corriente.append(I)
    plt.plot(voltaje,corriente)
    corriente=[]

plt.show()