%%% This script generates a look-up table for front wheel angles for given
%%% steering input angles 


% housekeeping
close all
clear
clc

% Steering angle
input_steering_angle = [-124:1:124]; % [deg] range of motion of pinion gear
rack_min = -41.275; % [mm] left most position of rack 
rack_max = 41.275; % [mm] right most position of rack
rack_pos = linspace(rack_min, rack_max, length(input_steering_angle)); % [mm] array of rack positions

% Input geometry
a = 355; % Tie rod length [mm]
b = 92.5; % knuckle to pivot length [mm]
x = 52; % [mm]
y_centre = 328; % [mm]

% Calculate wheel angle at each rack position value
for i = [1:1:length(rack_pos)]
     
% FRONT RIGHT WHEEL
    y_r(i) = y_centre-rack_pos(i);

    z_r(i) = sqrt(x^2 + y_r(i)^2); 
    phi_r(i) = atan(x/y_r(i)); 
    theta_r(i) = acos((z_r(i)^2 + b^2 - a^2)/(2*b*z_r(i))); % cosine rule
    wheel_angle_r(i) = theta_r(i) + phi_r(i); 
    wheel_angle_deg_r(i) = wheel_angle_r(i)*180/pi - 105.720479102835; % 105.72 deg is the 'resting' angle which gives zero wheel angle at zero steering angle

% FRONT LEFT WHEEL 
    wheel_angle_deg_l(i)= - wheel_angle_deg_r(i);

end

wheel_angle_deg_l = fliplr(wheel_angle_deg_l); % flip array from left to right due to steering assembly symmetry 

% Generate plot of look-up table for front left and right wheel angles
% against input steering angle
plot(input_steering_angle, wheel_angle_deg_l, 'b')
hold on
plot(input_steering_angle, wheel_angle_deg_r, 'r')
hold off
ylabel("Front wheel angle, deg");
xlabel("Input steering angle, deg");
legend('Left front wheel', 'Right front wheel', 'location', 'northwest')