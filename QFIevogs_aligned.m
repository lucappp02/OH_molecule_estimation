function y = QFIevogs_aligned(E,index,t)
h = 1.055e-34;
delta = 2*pi*1.667*10^9;
Kb = 1.38e-23;
me = 1.66 *3.336e-30;
E = me *E/Kb;
D = h*delta/Kb;
t = t/h*Kb;
o = sqrt(D^2+36/25*E^2);
x = sqrt(D^2+4/25*E^2);
Co = (cos(o*t)-1)/x^2;
Cx = (cos(x*t)-1)/x^2;
So = (sin(o*t)-o*t)/o^3;
Sx = (sin(x*t)-x*t)/x^3;
if index ==1
    y = (225* (5* t^2 + 8 *(2* Co^2 + So* t) *D^2 +  16* So^2* D^4) + 16* (Sx - 9* So)^2 *D^2 *(E)^2)/(1296* t^2* (16* Co^2* D^2 + (t + 4 *So* D^2)^2));
else
    y = (625* (25* t^2 + 4* (Sx + 9 *So) *t *D^2 + 8* D^2 *(Cx^2 + 9* Co^2 + (Sx^2 + 9* So^2)* D^2)) +  800* (Sx^2 + 81 *So^2)* D^2* (E)^2)/(16* t^2* (125* (5 *t^2 + 4* (Sx + 9* So) *t *D^2 +  8 *D^2 *(Cx^2 + 9* Co^2 + (Sx^2 + 9 *So^2)* D^2)) + 96* (Sx^2 + 12 *Sx* So + 81* So^2) *D^2 *(E)^2));
end
end