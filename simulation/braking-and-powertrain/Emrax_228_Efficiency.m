[Motor_Speed,Normalised_Motor_Speed,Motor_Torque,Normalised_Motor_Torque,Efficiency] = readvars('EMRAX_Efficiency');   

N = 25;
xvec = linspace(min(Motor_Speed), max(Motor_Speed), N);
yvec = linspace(min(Motor_Torque), max(Motor_Torque), N);
[X, Y] = ndgrid(xvec, yvec);
F = scatteredInterpolant(Motor_Speed, Motor_Torque, Efficiency);
Z = F(X, Y);
surf(X, Y, Z);
s.EdgeColor = 'none';

%colormap summer
title('EMRAX 228 Efficiency')
xlabel('Speed /rpm') 
ylabel ('Torque /Nm')
zlabel ('Efficiency')
