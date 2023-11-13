import pygad
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


#data = {'Ir': Ir, 'Pmax': Pmax, 'Vmax': Vmax, 'Vvec':Vvec, 'Pvec': P, 'Ivec': Ivec} lo que entrega simulate2
#ARRAYS ARE CREATED TO SAVE THE POWER, VOLTAGE AND CURRENT CURVES.
vvec=np.ones(309)
ivec=np.ones(309)
pvec=np.ones(309)
pmax=np.zeros(309)
iter=0
iter2=0
voc=44.2
irr =np.array([[1000.,1000.,1000.],[1000.,500.,200.],[100.,100.,100.]])#IRRADIANCE PROFILES USED.
T=25.0
vgbest=0
fitness=0


def fitness_func(solution,solution_idx):
    #OBJECTIVE FUNCTION, RECEIVES THE NEW SOLUTIONS BY COMPARING THEM WITH THE OBJECTIVE VALUE.
    global current_irr,fitness,iter,iter2,vvec,ivec,vgbest
    paramss = np.append(current_irr,[T,solution])
    dataa = pvsim(paramss.tolist())
    ipv=dataa['Ir']
    pvv=ipv*solution
    if(pvv[0]>fitness):
        vgbest=solution
    fitness = pvv[0]
    vvec[iter]=solution
    ivec[iter2]=ipv
    iter=iter+1
    iter2=iter2+1
    return fitness


def po(vact,vpas,ipas,step):
    #FUNCTION THAT PLAYS THE P&O
    global current_irr

    params = np.append(current_irr,[T,vact])
    data = pvsim(params.tolist())
    I_pv=data['Ir']
    pact=vact*I_pv
    ppas=vpas*ipas

    
    if pact > ppas:                               #SI P_ACT ES MAYOR A P_PAS
        if vact > vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #subir
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #bajar
            return vact,vpas,ipas
    else:                                            #SI P_ACT NO ES MAYOR A P_PAS
        if vact < vpas:
            vpas=vact
            ipas=I_pv
            vact=vact+step #subir
            return vact,vpas,ipas
        else:
            vpas=vact
            ipas=I_pv
            vact=vact-step #bajar
            return vact,vpas,ipas


for instance in range(3):
    num=0
    num2=0
    current_irr=irr[instance,:] #CURRENT IRRADIANCE
    p_creature=[[0.15*voc],[0.5*voc],[0.85*voc]] #INITIAL POSITION OF THE INDIVIDUALS

    #GA STAGE
    ga_instance = pygad.GA(num_generations=25,
                       num_parents_mating=2,
                       fitness_func=fitness_func,
                       initial_population=p_creature,
                       parent_selection_type='rank',
                       keep_parents=2,
                       crossover_type='single_point',
                       crossover_probability=0.8,
                       mutation_type='random',
                       mutation_num_genes=1,
                       random_mutation_min_val=-1,
                       random_mutation_max_val=1)

    ga_instance.run()

    step=0.2
    vpas=vvec[iter-2]
    ipas=1
    vact,vpas,ipas=po(vgbest,vpas,ipas,step) #INITIALIZE P&O

    vvec[iter]=vact
    ivec[iter2]=ipas
    iter=iter+1
    iter2=iter2+1

    #P&O STAGE
    for num in range(73):

        vact,vpas,ipas=po(vact,vpas,ipas,step)
        vvec[iter]=vact
        ivec[iter2]=ipas
        iter=iter+1
        iter2=iter2+1

x=np.arange(0,309*(6/305),6/305)

#FUNTION TO PLOT MAXIMUM POWER IN TIME.
for z in range(309):
    if(z<102):
        pmax[z]=380
    if(z>=102 and z<204):
        pmax[z]=135
    if(z>=204):
        pmax[z]=33

fig, axes = plt.subplots(nrows=3, ncols=1)


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

#plt.savefig("gapo.pdf")
plt.show()
