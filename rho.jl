function rho(B, E, theta, timevar, Temp)
   
    hbar = 1.055e-34
    delta = 2 * π * 1.667e9
    Kb = 1.38e-23
    timevar = timevar / hbar * Kb


    # Initial thermal state
    
    H = hamOH(B, E, theta)
    H0 = hamOH(0.,0.,1.);
    rho_0 = exp(-H0/Temp)/tr(exp(-H0/Temp))

    #unitary evolution
    U = exp(-1im * timevar * H)

    # Evolved state
    y = U * rho_0 * U'

    return y
end