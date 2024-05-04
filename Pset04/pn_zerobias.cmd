File {
    * Inputs
    Grid      = "pn_msh.tdr"

    * Outputs
    Plot      = "pn_zerobias.tdr"
    Current   = "pn_zerobias.plt"
    Output    = "pn_zerobias.log"
}

Electrode {
    { Name="anode"    Voltage=0.0 }
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
}
















