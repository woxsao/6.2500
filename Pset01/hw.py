import numpy as np


#CONSTANTS =============================================================================================
epsilon_0 = 8.85*10**(-12)
epsilon_SiO2 = 3.9
d = 5*10**(-9)
width = 200*10**(-9)
length = 2*10**(-6)
Cin = 2*epsilon_0*epsilon_SiO2*width*length/(d)
Cout = 1/4*Cin

Vdd = 1
Vt = 0.4
Ion = 2*0.0004*(Vdd-Vt)**2
Ron = Vdd/Ion


print("Cin", Cin)
print("Cout", Cout)

#QUESTION 1A=============================================================================================
print("QUESTION 1A ========================================================")
Rw = 1200
Cw = 0.15*10**(-15)*60
print("Cw", Cw)
FO = 2
LD = 18
tclk = LD*1/2*np.log(2)*(Ron*(Cout + Cw+2*Cin)+Rw*(Cw/2+FO*Cin))
print("tclk", tclk)
print("freq", 1/tclk)


#QUESTION 1B=============================================================================================
print("QUESTION 1B ========================================================")
#NUMBER OF GATES PER STAGE CALCULATION
gate_per_stage = 0
for i in range(LD):
    coeff = 2**i
    gate_per_stage += coeff

edgate = (Cout + Cw + 2*Cin)*Vdd**2*0.1
print("edgate", edgate)
edtotal = gate_per_stage*5*edgate
print("edtotal", edtotal)
Ioff = 2*1*np.exp(-0.4/0.025)
Roff = Vdd/Ioff

Eleakgate = Vdd**2/Roff*tclk
print("eleakgate", Eleakgate)
eleaktotal= gate_per_stage*5*Eleakgate
print("eleak total", eleaktotal)
etotal = eleaktotal + edtotal
print("energy total ", etotal)

print("percent leak", eleaktotal/etotal)
print("percent dynamic", edtotal/etotal)
edp = etotal*tclk
print("Edp", edp)


#QUESTION 2A =============================================================================================
print("QUESTION 2A ========================================================")
tclk_graphene = LD*1/2*np.log(2)*(Ron*(Cout + Cw/2+FO*Cin)+Rw/100*(Cw/4+FO*Cin))
print("tclk graphene: ", tclk_graphene)
print("frequency graphene: ", 1/tclk_graphene)
freq_graphene = 1/tclk_graphene
print(tclk-tclk_graphene)

edynamic_graphene = gate_per_stage*5*(Cout + Cw/2+2*Cin)*Vdd**2*0.1
eleak_graphene = gate_per_stage*5*Vdd**2/Roff*tclk_graphene
print(edynamic_graphene)
total_energy_graphene = edynamic_graphene + eleak_graphene
print("graphene total energy", total_energy_graphene)
print("improvement: ", (total_energy_graphene- etotal)/etotal)

edp_graphene = tclk_graphene*total_energy_graphene
print("EDP Graphene", edp_graphene)



#QUESTION 2B=============================================================================================
print("QUESTION 2B ========================================================")
Vt_cnt = 0.4
Vdd_cnt = 0.8

Ioff_cnt = 2*np.exp(-Vt_cnt/0.025)
Roff_cnt = Vdd_cnt/Ioff_cnt
Ion_cnt = 2*0.0004*(Vdd_cnt-Vt_cnt)**2
Ron_cnt = Vdd_cnt/(2*Ion_cnt)
Cin_cnt = 2*epsilon_0*epsilon_SiO2*(40*10**(-9)*length)/d
Cout_cnt = 1/4*Cin_cnt
tclk_cnt = LD * 1/2*np.log(2)*(Ron_cnt*(Cout_cnt + Cw+2*Cin_cnt)+Rw*(Cw/2+2*Cin_cnt))
print("tclk cnt:", tclk_cnt)
print("cnt freq:", 1/tclk_cnt)


etotal_cnt = 5*gate_per_stage*(Vdd_cnt**2*0.1*(Cout_cnt+Cw+FO*Cin_cnt)+Vdd_cnt**2/Roff_cnt*tclk_cnt)
print("etotal cnt", etotal_cnt)
print("cnt improvement:", (etotal_cnt-etotal)/etotal)
edp_cnt = etotal_cnt*tclk_cnt
print("EDP CNT", edp_cnt)


#Question 2D=============================================================================================
print("QUESTION 2D ========================================================")
minimize = True
min_edp = np.inf
min_vdd = np.inf
min_vt = np.inf
min_tclk = 0
v_dd_range = np.arange(0.3,1.0001,0.001)



if(minimize):
    for v_dd in v_dd_range:
        v_t_range = np.arange(0.1,v_dd-0.1,0.001)
        for v_t in v_t_range:
            cur_Ion = 2*0.0004*(v_dd-v_t)**2
            cur_Ron_cnt = v_dd/(2*cur_Ion)

            cur_Ioff = 2*np.exp(-v_t/0.025)
            cur_Roff_cnt = v_dd/(cur_Ioff)

            cur_tclk = LD*1/2*np.log(2)*(cur_Ron_cnt*(Cout_cnt + Cw+FO*Cin_cnt)+Rw*(Cw/2+FO*Cin_cnt))
            cur_etotal = 5*gate_per_stage*(v_dd**2*0.1*(Cout_cnt+Cw+FO*Cin_cnt)+v_dd**2/cur_Roff_cnt*cur_tclk)
            cur_edp = cur_tclk*cur_etotal
            
            if(cur_edp < min_edp):
                min_edp = cur_edp
                min_vdd = v_dd
                min_vt = v_t
                min_tclk = cur_tclk
    print("minimized edp:", min_edp)
    print("minimized vdd:", min_vdd)
    print("minimized vt:", min_vt)
    print("minimized tclk:", min_tclk)
    print("minimized freq:", 1/min_tclk)


#QUESTION 3A =========================================================================================
print("QUESTION 3A ========================================================")
Cin = 2*epsilon_0*epsilon_SiO2*width*length/(d)
Cout = 1/4*Cin

Vdd = 1
Vt = 0.4
Ion = 2*0.0004*(Vdd-Vt)**2
Ron = Vdd/Ion


Rw = 1200
Cw = 0.15*10**(-15)*60
print("Cw", Cw)
FO = 2
LD = 18
tclk = LD*1/2*np.log(2)*(Ron*(Cout + 100*Cw+2*Cin)+100*Rw*(Cw*50+FO*Cin))
print("tclk", tclk)
print("freq", 1/tclk)


edgate = (Cout + 100*Cw + 2*Cin)*Vdd**2*0.1
print("edgate", edgate)
edtotal = gate_per_stage*5*edgate
print("edtotal", edtotal)
Ioff = 2*1*np.exp(-0.4/0.025)
Roff = Vdd/Ioff

Eleakgate = Vdd**2/Roff*tclk
print("eleakgate", Eleakgate)
eleaktotal= gate_per_stage*5*Eleakgate
print("eleak total", eleaktotal)
etotal = eleaktotal + edtotal
print("energy total ", etotal)

print("percent leak", eleaktotal/etotal)
print("percent dynamic", edtotal/etotal)
edp = etotal*tclk
print("Edp", edp)


#QUESTION 3B =============================================================================================
print("QUESTION 3B ========================================================")
minimize = True
min_edp = np.inf
min_vdd = np.inf
min_vt = np.inf
min_tclk = 0
v_dd_range = np.arange(0.3,1.0001,0.001)
min_energy = 0
if(minimize):
    for v_dd in v_dd_range:
        v_t_range = np.arange(0.1,v_dd-0.1,0.001)
        for v_t in v_t_range:
            cur_Ion = 2*0.0004*(v_dd-v_t)**2
            cur_Ron_cnt = v_dd/(2*cur_Ion)

            cur_Ioff = 2*np.exp(-v_t/0.025)
            cur_Roff_cnt = v_dd/(cur_Ioff)

            cur_tclk = LD*1/2*np.log(2)*(cur_Ron_cnt*(Cout_cnt + 100*Cw+FO*Cin_cnt)+100*Rw*(100*Cw/2+FO*Cin_cnt))
            cur_etotal = 5*gate_per_stage*(v_dd**2*0.1*(Cout_cnt+100*Cw+FO*Cin_cnt)+v_dd**2/cur_Roff_cnt*cur_tclk)
            cur_edp = cur_tclk*cur_etotal
            
            if(cur_edp < min_edp):
                min_edp = cur_edp
                min_vdd = v_dd
                min_vt = v_t
                min_tclk = cur_tclk
                min_energy = cur_etotal
    print("minimized edp:", min_edp)
    print("minimized vdd:", min_vdd)
    print("minimized vt:", min_vt)
    print("minimized tclk:", min_tclk)
    print('minimized freq:', 1/min_tclk)
    print("minimized energy:", min_energy)
