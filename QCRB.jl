function QCRB(B, E, theta, timevar, Temp)
    hbar = 1.055e-34
    delta = 2 * π * 1.667e9
    Kb = 1.38e-23
    timevar = timevar / hbar * Kb

    eta = (hbar * delta / (2 * Kb * Temp))
    rho_0 = exp(eta .* Diagonal([1, 1, 1, 1, -1, -1, -1, -1])) ./ (8 * cosh(eta))

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

    # Compute Hop as a 3D array
    Hop  = Matrix{ComplexF64}[]
    U = exp(-1im * timevar * H)
    
    for m in 1:3
        
        push!(Hop, -1im *U'* dU[m])
        
    end
    
    Q = Matrix{ComplexF64}
    
    Q = zeros(ComplexF64, 3, 3)

    for m in 1:3
        for n in 1:3
            for i in 1:4
                for j in 5:8
                   
                    Q[m, n] += real((Hop[m][i, j]) .* (Hop[n][j, i]))
                   
                end
            end
        end
    end
    Q = real(Q)




    y = (1/(tanh(eta)^2))*tr(inv(Q))
    return y
end