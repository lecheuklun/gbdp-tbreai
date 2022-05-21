% CONVERTDATA converts TTC data to format compatible for MFEval solver. 

clc; clear;

load B1965run24.mat;

% load('/Users/16koc1/Desktop/matlab/MFeval/doc/examples/TyreData.mat', 'TyreData')

m = length(ET);

RR = ones(m,1) * 16 ./ 39.37; % Rolling radius (m)

Time = ET; % 1 (s)
Vx = V ./ 3.6; % 2 Long Speed (m/s)
W = Vx ./ RR; % 3 Rotational speed (rad/s)
SA = SA .* (pi/30); % 4 Slip angle (rad)
IA = IA .* (pi/30); % 5 Inclination angle (rad)
Rl = RR; % 6 Loaded radius (m)
Re = RR; % 7 Effective rolling radius (m)
P = P * 1000; % 8 Pressure (Pa)
Fx = FX; % 9
Fy = FY; % 10
Fz = -FZ; % 11 Normal force (N)
Mx = MX;
Mz = MZ;
SR = SL;
Phit = zeros(m,1);

TyreData = table(Time,Vx,W,SA,IA,Rl,Re,P,Fx,Fy,Fz,Mx,Mz,SR,Phit);

save('TyreDataRear.mat','TyreData')
