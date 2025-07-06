function tight_bound(B, E, theta, timevar,Temp)
    hbar = 1.055e-34;
    Delta = 2*pi*1.667*10^9*hbar;
    Kb = 1.38e-23;
    timevar = timevar/hbar*Kb;
    prefactor = tanh(Delta/(2*Temp*Kb))^2;

    # Define matrices DBH, DExH, DEyH
    DBH = Diagonal([-6/5, -2/5, 2/5, 6/5, -6/5, -2/5, 2/5, 6/5])
    
    DEzH = [
        0  0  0  0  3/5  0  0  0;
        0  0  0  0  0  1/5  0  0;
        0  0  0  0  0  0  -1/5  0;
        0  0  0  0  0  0  0  -3/5;
        3/5  0  0  0  0  0  0  0;
        0  1/5  0  0  0  0  0  0;
        0  0  -1/5  0  0  0  0  0;
        0  0  0  -3/5  0  0  0  0
    ]

    DEyH = [
        0  0  0  0  0  -sqrt(3)/5  0  0;
        0  0  0  0  -sqrt(3)/5  0  -2/5  0;
        0  0  0  0  0  -2/5  0  -sqrt(3)/5;
        0  0  0  0  0  0  -sqrt(3)/5  0;
        0  -sqrt(3)/5  0  0  0  0  0  0;
        -sqrt(3)/5  0  -2/5  0  0  0  0  0;
        0  -2/5  0  -sqrt(3)/5  0  0  0  0;
        0  0  -sqrt(3)/5  0  0  0  0  0
    ]
    DH = Matrix{ComplexF64}[]
    push!(DH,DBH)
    push!(DH,DEzH)
    push!(DH,DEyH)

    H = hamOH(B, E, theta)
    E, en_eig = eigen(H)

    

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
    dummy = zeros(ComplexF64, 8,8)
    for m in 1:3
        dummy = 1im * dU[m]' * U
        push!(Hop, dummy)
        dummy = zeros(ComplexF64, 8,8)
    end

    # Initialize U and Q
    U_res = Matrix{Float64}
    Q = Matrix{Float64}
    U_res = zeros(Complex{Float64}, 3, 3)
    Q = zeros(Complex{Float64}, 3, 3)
    dummy1 = Matrix{ComplexF64}
    dummy1 = zeros(ComplexF64, 8,8)
    dummy = zeros(ComplexF64, 8,8)
    for m in 1:3
        for n in 1:3
            for i in 1:4
                for j in 5:8
                    dummy = Hop[m]
                    dummy1 = Hop[n]

                    U_res[m, n] += imag(dummy[i, j] * dummy1[j, i])
                    Q[m, n] += real(dummy[i, j] * dummy1[j, i])
                    dummy1 = zeros(ComplexF64, 8,8)
                    dummy = zeros(ComplexF64, 8,8)
                end
            end
        end
    end

    # Return the maximum eigenvalue of Q^(-1) * U_res
    operator = inv(Q)*U_res*inv(Q)
    operator_squared = operator'*operator
    g,av = eigen(operator_squared)
    y = sum(g.^(1/2)./prefactor)


    return y
end
