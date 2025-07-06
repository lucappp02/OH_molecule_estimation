function y = hamOH(B, E,t)
h = 1.055e-34;
delta = 2*pi*1.667*10^9;
Kb = 1.38e-23;

me = 1.66 *3.336e-30;
mb = 9.274e-24;
B= mb*B/Kb;
E = me *E/Kb;
D = h*delta/Kb;
y = [-1/2*D-6/5*B,0,0,0, 3/5*E*cos(t), -sqrt(3)/5*E*sin(t), 0,0;
    0, -1/2*D-2/5*B,0,0, -sqrt(3)/5*E*sin(t), 1/5*E*cos(t), -2/5*E*sin(t), 0;
    0,0,-1/2*D+2/5*B, 0 ,0, -2/5*E*sin(t), -1/5*E*cos(t), -sqrt(3)/5*E*sin(t);
    0,0,0,-1/2*D+6/5*B, 0,0, -sqrt(3)/5*E*sin(t), -3/5*E*cos(t);
    3/5*E*cos(t), -sqrt(3)/5*E*sin(t), 0,0,1/2*D-6/5*B,0,0,0;
    -sqrt(3)/5*E*sin(t), 1/5*E*cos(t),-2/5*E*sin(t), 0, 0, 1/2*D-2/5*B, 0, 0;
    0, -2/5*E*sin(t), -1/5*E*cos(t), -sqrt(3)/5*E*sin(t), 0, 0, 1/2*D+2/5*B, 0;
    0, 0, -sqrt(3)/5*E*sin(t), -3/5*E*cos(t), 0, 0,0, 1/2*D+6/5*B];
end
    