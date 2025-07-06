function y = R(B, E, t,time)



h = 1.055e-34;

kb = 1.38e-23;
time = time/h*kb;
DBH=  diag([-6/5, -2/5, 2/5, 6/5, -6/5, -2/5, 2/5, 6/5]);

DExH= [0,0,0,0, 3/5, 0, 0,0;
    0, 0,0,0, 0, 1/5, 0, 0;
    0,0,0, 0 ,0, 0, -1/5, 0;
    0,0,0,0, 0,0,0, -3/5;
    3/5, 0, 0,0,0,0,0,0;
   0, 1/5,0, 0, 0,0, 0, 0;
    0, 0, -1/5, 0, 0, 0,0, 0;
    0, 0, 0, -3/5, 0, 0,0, 0];


 DEyH=  [0,0,0,0,0, -sqrt(3)/5, 0,0;
    0, 0,0,0, -sqrt(3)/5, 0, -2/5, 0;
    0,0,0, 0 ,0, -2/5, 0, -sqrt(3)/5;
    0,0,0,0, 0,0, -sqrt(3)/5, 0;
    0, -sqrt(3)/5, 0,0,0,0,0,0;
    -sqrt(3)/5, 0,-2/5, 0, 0,0, 0, 0;
    0, -2/5, 0, -sqrt(3)/5, 0, 0,0, 0;
    0, 0, -sqrt(3)/5, 0, 0, 0,0, 0];

DH(:,:,1) = DBH;
DH(:,:,2) = DExH;
DH(:,:,3) = DEyH;


H = hamOH(B,E,t);
[V,D] = eig(H);
d = diag(D);
% U = zeros(8,8);
% for i = 1:8
%     for j = 1:8
%         for k = 1:8
%             U(i,j) = U(i,j)+ exp(-sqrt(-1)*time*d(k))*V(i,k)*V(j,k);
%         end
%     end
% end
dU = zeros(8,8,3);
for n = 1:3
    for i = 1:8
        for j = 1:8
            for l = 1:8
                for m = 1:8
                    if l~= m
                        dU(i,j,n) = dU(i,j,n) -2* sqrt(-1)*V(i,m)*V(j,l).*V(:,m)'*DH(:,:,n)*V(:,l).*exp(-sqrt(-1)*time*(d(l)+d(m))/2)*sin(time*(d(l)-d(m))/2)/(d(l)-d(m)); %CONTROLLARE SEGNI
                    else
                       dU(i,j,n) = dU(i,j,n) - sqrt(-1)*time.*V(:,m)'*DH(:,:,n)*V(:,m).*V(j,m).*V(i,m).*exp(-sqrt(-1)*time*d(m));
                    end
                end
            end
        end
    end
end
Hop = zeros(8,8,3);
U = expm(-1i*time*H );
for m = 1:3
    Hop(:,:,m) = 1i*dU(:,:,m)'*U;
end






U = zeros(3,3);
for m = 1:3
    for n = 1:3
        for i = 1: 4
            for j = 5:8
                U(m,n) = U(m,n)+imag(Hop(i,j,m)*Hop(j,i,n));
            end
        end
    end
end

Q= zeros(3,3);

for m = 1:3
    for n = 1:3
        for i = 1:4
            for j = 5:8
                Q(m,n) = Q(m,n)+real( Hop(i,j,m)*Hop(j,i,n));
            end
        end
    end
end

y = max(abs(eig(1i*Q^(-1)*U)));
end