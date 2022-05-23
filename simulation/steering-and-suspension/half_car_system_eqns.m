function [dz, a1, a2] = half_car_system_eqns(t,z)
% Half car model
% This function calculates the system equations using a state space model
% state vector Z = [z1;z2;z3;z4;z5;z6;z7;z8] = [x;theta;x1;x2;x_dot;theta_dot;x1_dot;x2_dot]
% state derivative vector dZ = [x_dot;theta_dot;x1_dot;x2_dot;x_dot_dot;theta_dot_dot;x1_dot_dot;x2_dot_dot]

% System Parameters
% (measured)
M   = 207/2;         % Half of body mass                             [kg]   
a1  = 0.8415;        % Distance of CoG from front axle               [m]
a2  = 0.6885;        % Distance of CoG from rear axle                [m]

% (calculated using lumped mass model)
m1  = M*a2/(a1+a2);  % Mass of a front unsprung mass                 [kg]
m2  = M*a1/(a1+a2);  % Mass of a rear unsprung mass                  [kg] 
I  = M*a1*a2;        % Half of body lateral mass moment of inertia   [kg.m^2]

% (unknown for car, so values obtained from research of similar formula student cars to give best approximation possible)
k1  = 4480;          % Spring constant suspension front              [N/m]
k2  = 5355;          % Spring constant suspension rear               [N/m]
c1  = 260;           % Damping constant suspension front             [N.s/m]
c2  = 390;           % Damping constant suspension rear              [N.s/m]
kt1 = 125000;        % Spring constant tire front                    [N/m]
kt2 = 125000;        % Spring constant tire rear                     [N/m]

% System Matrix
A = [0 0 0 0 1 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1;
     -(k1+k2)/M        (k1*a1-k2*a2)/M           k1/M           k2/M           -(c1+c2)/M        (c1*a1-c2*a2)/M            c1/M        c2/M    ;
     (a1*k1-a2*k2)/I   -((a1^2)*k1+(a2^2)*k2)/I  -a1*k1/I       a2*k2/I        (a1*c1-a2*c2)/I   -((a1^2)*c1+(a2^2)*c2)/I   -a1*c1/I    a2*c2/I ;
     k1/m1             -k1*a1/m1                 -(k1+kt1)/m1   0              c1/m1             -c1*a1/m1                  -c1/m1      0       ;
     k2/m2             k2*a2/m2                  0              -(k2+kt2)/m2   c2/m2             c2*a2/m2                   0           -c2/m2  ];

% State Derivative = System Matrix * State Vector
dz = A * [z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8)];