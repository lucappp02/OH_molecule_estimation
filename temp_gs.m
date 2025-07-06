%% intervallo di confidenza per il GS
E_min = 1e5;
E_max = 1e6;
Bmin = 1;
Bmax = 1;
theta_min = 0;
theta_max = pi/2;
mesh = 10^2;
step_T = 1e-4;
T0 = 1e-4;
step_theta = (theta_max-theta_min)/mesh;
step_B = (Bmax-Bmin)/mesh;
step_E = (E_max-E_min)/mesh;
confidence = 0.9;
dummy = 1000;

    B = Bmin
    for j = 1:mesh
        E = E_min+j*step_E;
        for k=0:mesh-1
            theta = theta_min+k*step_theta;
            
            %d= sort(eig(hamOH(B, E, theta)));
            
            T = T0;
            l = 0;
            P = 1;
            while (P >=confidence || isnan(P))&&T<dummy
                    
                T = T+l*step_T;
               d =  sort(eig(expm(-hamOH(B,E,theta)/T)/trace(expm(-hamOH(B,E,theta)/T))));
               P = d(8);
                %P = exp(-2*d(1)/T)/(sum(exp(-d/T)))^2;
            l = l+1;
            end
            
            
            if T <dummy
                dummy = T;
               
                save_j = j;
                save_k = k;
                save_P = P
            end
        end
        j
    end
    

dummy
 B =  Bmin
  E = E_min+save_j*step_E
  t = theta_min+save_k*step_theta
 d =  sort(eig(expm(-hamOH(B,E,t)/dummy)/trace(expm(-hamOH(B,E,t)/dummy))));
 check = save_P-d(8)^2