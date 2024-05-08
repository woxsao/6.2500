File { 
* The File section defines the input and output files of the simulation
* Input Files
Grid = "Si_MOSFET_msh.tdr"   *fix the FILENAME

* Output Files
Current = "FILENAME_def.plt" 
Plot = "FILENAME_def.tdr" 
Output = "FILENAME_def.log" 
}

Electrode {
{ Name="Drain_contact" Voltage=2.7 } * CHANGE THIS TO SET VDD
{ Name="Source_contact" Voltage=0.0 } 
{ Name="Gate_contact" Voltage=0.0 }
{ Name="Body_contact" Voltage=0.0 }
}

Physics { 
Mobility( DopingDep HighFieldSat Enormal )
EffectiveIntrinsicDensity ( OldSlotBoom )
} Plot { eDensity hDensity eCurrent hCurrent Potential SpaceCharge ElectricField eMobility hMobility eVelocity hVelocity
Doping DonorConcentration AcceptorConcentration }

Math {
Extrapolate 
RelErrControl *on by default
Iterations=50 
Notdamped=100 
}

Solve {
Coupled(Iterations=100){ Poisson }
Coupled{ Poisson Electron Hole }
*-Bias Cathode to target bias
Quasistationary(
InitialStep=0.001 Increment=1.1
MinStep=1e-5 MaxStep=0.05
Goal{ Name="Gate_contact" Voltage=2.7 } *CHANGE THIS TO SET VDD
){ Coupled{ Poisson Electron Hole }}
} 