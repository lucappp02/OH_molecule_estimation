using QuanEstimation
using LinearAlgebra
using Plots
using LaTeXStrings

mesh = 101;
tmin = 1e-10;
tmax = 1e-9;
x = range(tmin, tmax, mesh)




Eval = 1e5;
Bval = 0.1;
Tempval = 1;
thetaval = pi/8;



y3 = zeros(mesh)
z3 = zeros(mesh)
w3 = zeros(mesh)

for j=1:mesh
    local_rho_mat  = rho(Bval, Eval, thetaval,x[j],Tempval);
    local_drho_mats = derivative_rho(Bval, Eval, thetaval,x[j],Tempval);
    local_Wmat = [
            1	0	0;
            0	1	0	;
            0	0	1	
            ]
    local_hcrbval = HCRB(local_rho_mat,local_drho_mats,local_Wmat)
    qcrbval = QCRB(Bval, Eval, thetaval,x[j],Tempval)

    y3[j] = local_hcrbval/qcrbval-1;
    z3[j] = R(Bval, Eval, thetaval,x[j])
    w3[j] = tight_bound(Bval, Eval, thetaval,x[j],Tempval)/qcrbval
    
    
   

end
p3 = plot(x,y3,label=L"\frac{C^H}{C^S}-1",title= L" B=0.1\, T\;, E=10^5\, kV/m\;, \theta=\frac{\pi}{8}\;, T=1\, K", xlabel= L"t")
plot!(x,z3, label = L"R")
plot!(x,w3, label = L"‖Q^{-1}DQ^{-1}‖_1/C^S")

Eval = 1e5;
Bval = 0.1;
Tempval = 1;
thetaval = pi/3;



y4 = zeros(mesh)
z4 = zeros(mesh)
w4 = zeros(mesh)
for j=1:mesh
    local_rho_mat  = rho(Bval, Eval, thetaval,x[j],Tempval);
    local_drho_mats = derivative_rho(Bval, Eval, thetaval,x[j],Tempval);
    local_Wmat = [
            1	0	0;
            0	1	0	;
            0	0	1	
            ]
    local_hcrbval = HCRB(local_rho_mat,local_drho_mats,local_Wmat)
    qcrbval = QCRB(Bval, Eval, thetaval,x[j],Tempval)

    z4[j] = R(Bval, Eval, thetaval,x[j])
   
    y4[j] = local_hcrbval/qcrbval-1;
    w4[j] = tight_bound(Bval, Eval, thetaval,x[j],Tempval)/qcrbval

end
p4 = plot(x,y4,label=L"\frac{C^H}{C^S}-1",title= L" B=0.1\, T \;,E=10^5\, kV/m\;, \theta=\frac{\pi}{3}\;, T=1\, K", xlabel= L"t")
plot!(x,z4, label = L"R")
plot!(x,w4, label = L"‖Q^{-1}DQ^{-1}‖_1/C^S")

plot(p3, p4,lw = 4 ,layout=(1,2),size=(1000,600))
