function [delta_x_front, delta_x_rear] = half_car_using_ode45(z0)
% This function outputs vertical displacements of the car body at the front
% and rear axles using the system of equations from the function "half_car_system_eqns.m"
% and MATLAB's ode45 function.
% These outputs are then called by the "kinematic_susp_model_2D.m" file 
% where they are used to simulate the suspesnion's behaviour and calculate camber angles.

odefun = @half_car_system_eqns;
t0 = 0; % Initial time (s)
tend = 3; % End time (s)
tspan = [t0, tend]; % Time array for simulation (s)

% z0 contains the initial conditions of the state vector in [m; rad; m; m; m/s; rad/s; m/s; m/s]
[t, z] = ode45(odefun,tspan,z0); % Use MATLAB's "ode45" function to solve the system of ODEs using RK 4,5 methods

[~ , a1, a2] = half_car_system_eqns(t,z);
% Calculating car body vertical displacements at the front and rear axles
% using:
% output state variable z1 (vertical displacement of car body's CoG)
% output state variable z2 (pitch angle of car body)
% a1 and a2 (distance from vehicle CoG and front and rear axles respectively)
delta_x_front = z(:,1) - a1*sin(z(:,2)); % m
delta_x_rear = z(:,1) + a2*sin(z(:,2)); % m

% Plot the results
figure(1)
% M --- car
plot(t,z(:,1),'b') % Vertical Displacement body, m
hold on
plot(t,z(:,2),'r') % Pitch Angle body, deg
hold on
% m1 --- front wheel
plot(t,z(:,3),'k') % Vertical Displacement front wheel, m
hold on
% m2 --- rear wheel
plot(t,z(:,4),'m') % Vertical Displacement rear wheel, m
hold off

xlim([t0 tend])
xlabel('Time, s')
ylabel('State Variables')
legend('Vertical Displacement body, m','Pitch Angle body, deg', 'Vertical Displacement front wheel, m','Vertical Displacement rear wheel, m')

figure(2)
% Vertical Displacement at front axle
plot(t,delta_x_front,'k')
hold on
% Vertical Displacement at rear axle
plot(t,delta_x_rear,'m')
hold off

xlim([t0 tend])
xlabel('Time, s')
ylabel('Vertical Displacement, m')
legend('Vertical Displacement at front axle, m','Vertical Displacement at rear axle, m')