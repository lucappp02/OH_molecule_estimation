function y = QFIExEygibbs(B, E,t, T)

%COERENTE CON  QFIgsxy




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



[V,D] = eig(hamOH(B,E,t));
d = diag(D);




b =1/T;
Z= sum(exp(-b*d));
dE = zeros(3,8);
for m = 1:3
   dE(m,:) = diag(V'*DH(:,:,m)*V);
end

dZ = -b*dE*(exp(-b*d));  %vettore 3x1

Q = zeros(3,3);
 

for m = 1:3
    for n = 1:3
        for i = 1:8
             for j = 1:8
                if i~=j
                     Q(m,n) = Q(m,n)+ exp(-b*d(i))/Z*(2*tanh(b*(d(i)-d(j))/2)/(d(i)-d(j)))^2*V(:,i)'*DH(:,:,m)*V(:,j)*V(:,j)'*DH(:,:,n)*V(:,i);
                end 
            end
        end
       for i = 1:8
            Q(m,n) = Q(m,n)+ exp(-b*d(i))/Z*  (b*dE(m,i) + dZ(m)/Z).*  (b*dE(n,i) + dZ(n)/Z);
            
        end
    end
end



y = trace(Q^(-1));

end