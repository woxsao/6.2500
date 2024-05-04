File {
    * Inputs
    Grid      = "pn_msh.tdr"

    * Outputs
    Plot      = "pn_iv.tdr"
    Current   = "pn_iv.plt"
    Output    = "pn_iv.log"
}

Electrode {
    { Name="anode"    Voltage=-2.0 }
    { Name="cathode"  Voltage=0.0 }
}

Physics {
    EffectiveIntrinsicDensity( OldSlotboom )
    Mobility( DopingDep HighFieldSat Enormal )
    Recombination( SRH( DopingDep ) )
}

Plot {
    eDensity hDensity
    TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
    eMobility hMobility
    eVelocity hVelocity
    eQuasiFermi hQuasiFermi
    eQuasiFermiPotential hQuasiFermiPotential

    ElectricField/Vector Potential SpaceCharge

    Doping DonorConcentration AcceptorConcentration

    SRH

    eGradQuasiFermi/Vector hGradQuasiFermi/Vector
    eEparallel hEparallel eENormal hENormal

    BandGap 
    BandGapNarrowing
    Affinity
    ConductionBand ValenceBand
    eQuantumPotential
}

Math {
    Extrapolate 
    RelErrControl *on by default
}

Solve {
    # initial solution, starting voltage
    Poisson
    Coupled(Iterations=100) { Poisson }
    Coupled { Poisson Electron Hole }
    
    # ramp anode and save solutions :
    Quasistationary(
        InitialStep=0.001 Increment=1.1
        MinStep=0.0001 MaxStep=0.01
        Goal{ Name="anode" voltage=2.0 } *CHANGE THIS TO SET VOLTAGE SWEEP TARGET
    ) { Coupled { Poisson Electron Hole } }
}
















