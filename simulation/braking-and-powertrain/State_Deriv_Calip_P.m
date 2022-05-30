function dP = State_Deriv_Calip_P(t,P_Calipers,X,P_Ref)

V = 0.1;      % ASSUMED VALUE-----Volume of caliper side
B = 0.9*10^4; % ASSUMED VALUE-----Research Further
Cd = 0.62;    % ASSUMED VALUE-----Research Further /-
d = 0.005;    % ASSUMED VALUE-----Area of annular orifice
%X = 0.001;   % Spool displacement from centre (Zero-lapped valve)
rho = 1070;   % Density of DOT 4 brake fluid

if X > 0
    dP = (B/V)*Cd*pi*d*X*((2*(P_Ref-P_Calipers))/rho)^(1/2);
else
    dP = (B/V)*Cd*pi*d*X*((2*(P_Calipers-P_Ref))/rho)^(1/2);
end
