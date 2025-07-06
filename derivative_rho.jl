function derivative_rho(B, E, theta, timevar, Temp)
    hbar = 1.055e-34
    delta = 2 * π * 1.667e9
    Kb = 1.38e-23
    timevar = timevar / hbar * Kb

    eta = (hbar * delta / (2 * Kb * Temp))

    H = hamOH(B, E, theta)
    H0 = hamOH(0.,0.,1.);
    rho_0 = exp(-H0/Temp)/tr(exp(-H0/Temp))

    # Hamiltonian derivatives
    DBH = Diagonal([-6/5, -2/5, 2/5, 6/5, -6/5, -2/5, 2/5, 6/5])
    
    DEzH = [0 0 0 0 3/5 0 0 0;
            0 0 0 0 0 1/5 0 0;
            0 0 0 0 0 0 -1/5 0;
            0 0 0 0 0 0 0 -3/5;
            3/5 0 0 0 0 0 0 0;
            0 1/5 0 0 0 0 0 0;
            0 0 -1/5 0 0 0 0 0;
            0 0 0 -3/5 0 0 0 0]
    
    DEyH = [0 0 0 0 0 -sqrt(3)/5 0 0;
            0 0 0 0 -sqrt(3)/5 0 -2/5 0;
            0 0 0 0 0 -2/5 0 -sqrt(3)/5;
            0 0 0 0 0 0 -sqrt(3)/5 0;
            0 -sqrt(3)/5 0 0 0 0 0 0;
            -sqrt(3)/5 0 -2/5 0 0 0 0 0;
            0 -2/5 0 -sqrt(3)/5 0 0 0 0;
            0 0 -sqrt(3)/5 0 0 0 0 0]
    
    DH = Matrix{ComplexF64}[]
    push!(DH,DBH)
    push!(DH,DEzH)
    push!(DH,DEyH)

    H = hamOH(B, E, theta)
    E, en_eig = eigen(H)


    # Derivatives of the Unitary operator
    dU = Matrix{ComplexF64}[]
    dummy = Matrix{ComplexF64}
    dummy = zeros(ComplexF64, 8,8)

    for n = 1:3
        for i = 1:8
            for j = 1:8
                for k = 1:8
                    for l = 1:8
                        if l == k
                            dummy[i, j] -= 1im*timevar*exp(-1im*timevar*E[k])*  (en_eig[:,k]'*DH[n]*en_eig[:,k]) * en_eig[i,k]*conj(en_eig[j,k])
                        else
                            dummy[i, j] += (en_eig[:,l]'*DH[n]*en_eig[:,k])  /  (E[l]-E[k])  *  (exp(-1im*timevar*E[l])-exp(-1im*timevar*E[k])) *  en_eig[i,l]*conj(en_eig[j,k])
                        end
                    end
                end
            end
        end
        push!(dU, dummy)
        dummy = zeros(ComplexF64, 8,8)
    end

    U = exp(-1im * timevar * H)
    dummy = zeros(ComplexF64, 8,8)
    # Derivative of the evolved state
    der_rho  = Matrix{ComplexF64}[]
    for k = 1:3
        dummy = dU[k] * rho_0 * U'
        dummy = dummy+dummy'
        push!(der_rho, dummy)
        dummy = zeros(ComplexF64, 8,8)
    end

    return der_rho
end