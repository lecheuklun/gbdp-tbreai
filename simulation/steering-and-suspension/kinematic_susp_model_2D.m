%% This script models the front and rear suspension of the TBReAI FS car 
% using a 2D kinematic analysis.
%
% The script requires 2 user inputs:
% 1. Whether you would like to analyse the front or rear suspension
% 2. Whether you would like to run the simulation for the suspension's full
%    range of motion, or whether you would like to perform the simulation
%    using inputs from a half car model for which the initial conditions can
%    be controlled.

%% house keeping
close all
clear
clc

%% Get user input
% Front or rear suspension simulation
prompt = {'Would you like to simulate the front or rear suspension?'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'Enter "front" or "rear"'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

prompt2 = {'Would you like to simulate the full range of suspension motion, or the real-world scenario?'};
dlgtitle2 = 'Input';
dims2 = [1 50];
definput2 = {'Enter "full range" or "real-world"'};
answer2 = inputdlg(prompt2,dlgtitle2,dims2,definput2);

% Get vertical displacements of body at front and rear axles from half car model. 
if answer2 == "real-world"
    if answer == "front" 
        z0 = [0.2; 5*pi/180; 0; 0; 0; 0; 0; 0]; % arbitrary initial conditions in [m; rad; m; m; m/s; rad/s; m/s; m/s]
    else
       z0 = [0.1; 5*pi/180; 0; 0; 0; 0; 0; 0]; % arbitrary initial conditions in [m; rad; m; m; m/s; rad/s; m/s; m/s]
    end
       
    [delta_x_front, delta_x_rear] = half_car_using_ode45(z0); 
end

%% Inputs 
if answer == "front"
    % Front suspension data
    % Input lengths
    c = 310.6; % Upper link length (mm)
    a = 423.2; % Lower link length (mm)
    b = 141.2; % King pin length (mm)
    d = 110.6; % Chassis hinge plane vertical distance (mm)
    e = 4.0; % stub axle length (mm)
    delta = 106.5; % (mm)
    e2 = 115.0 + e; % length postion 2 (mm)
    theta2_design = 87.7*pi/180; % (rad) this is the value of theta_2 in the static resting position

    % Input angle
    alpha = 92.0*(pi/180); % King pin to SA angle (fixed) (rad)
    
    if answer2 == "real-world"
    % FROM HALF CAR MODEL
    % Vertical Displacement of car body at front axle (m)
    delta_x = delta_x_front;
    x_design = a*cos(theta2_design);
    theta2 = acos((x_design-delta_x*1000)/a); % have to convert delta_x_front from (m) to (mm)

    else 
    % FULL RANGE OF MOTION FOR ANIMATION
    % Motion of suspension
    theta2_bump1 = [87:-1:75]*(pi/180);
    theta2_rebound = [75:1:99]*(pi/180);
    theta2_bump2 = [99:-1:87]*(pi/180);
    theta2 = horzcat(theta2_bump1, theta2_rebound,theta2_bump2); % Lower link to chassis pivot plane angle (rad)
    end 

    % Front Wheel/Tyre Dimensions
    tyre_width = 190; % mm
    tyre_height = 406.4; % mm (same as diameter)

    
else
    % Rear suspension data
    % Input lengths
    c = 219.7; % Upper link length (mm)
    a = 285.2; % Lower link length (mm)
    b = 139.9; % King pin length (mm)
    d = 104.0; % Chassis hinge plane vertical distance (mm)
    e = 5.8; % stub axle length (mm)
    delta = 61.0; % (mm)
    e2 = 127.7 + e; % length postion 2 (mm)
    theta2_design = 90.5*pi/180; % (rad) this is the value of theta_2 in the static resting position
    
    % Input angle
    alpha = 95.1*(pi/180); % King pin to SA angle (fixed) (rad)

    if answer2 == "real-world"
    % FROM HALF CAR MODEL
    % Vertical Displacement of car body at rear axle (m)
    delta_x = delta_x_rear;
    x_design = a*cos(theta2_design);
    theta2 = acos((x_design-delta_x*1000)/a); % have to convert delta_x_rear from (m) to (mm)

    else 
    % FULL RANGE OF MOTION FOR ANIMATION
    % Motion of suspension
    theta2_bump1 = [90:-1:70]*(pi/180); 
    theta2_rebound = [70:1:110]*(pi/180);
    theta2_bump2 = [110:-1:90]*(pi/180);
    theta2 = horzcat(theta2_bump1, theta2_rebound, theta2_bump2); % Lower link to chassis pivot plane angle (rad)
    end 

    % Rear Wheel/Tyre Dimensions
    tyre_width = 208; % mm
    tyre_height = 406.4; % mm (same as diameter)

end 

% Apparent lengths
AN = sqrt( (a-delta)^2 + (d-delta.*cos(theta2)).^2 - 2*(a-delta)*(d-delta.*cos(theta2)).*cos(theta2) ); 

%% Outputs

% Calculated angles
beta = atan((sqrt(4*b^2*AN.^2 - (b^2 + AN.^2 - c^2).^2))./(b^2 + AN.^2 - c^2)); % AN to king pin angle (rad)
gamma = atan(((a-delta)*sin(theta2))./(d-a*cos(theta2))); % AN to vertical reference angle (rad)
phi = beta - gamma + alpha; % Stub axle to vertical refrence angle (rad)
camber = phi - pi/2; % Stub axle vertical face to reference plane angle (rad)

% Positions of linkage joints
A_x = a*sin(theta2); % Lower link outer
A_y = a*cos(theta2); 

B_x = A_x - b*sin(gamma-beta); % Upper link outer
B_y = A_y + b*cos(gamma-beta);

M_x = 0*theta2; % Lower link chassis
M_y = 0*theta2;

N_x = delta*sin(theta2); % Upper link chassis
N_y = (0 + d)*ones(1,length(theta2));

% Position of stub axle reference location
C_x = A_x + e*sin(phi) - (b/2)*sin(gamma-beta); % x plane position
C_y = A_y + e*cos(phi) + (b/2)*cos(gamma-beta); % y plane position
 
% Position of stub axle reference location 2
C_x2 = A_x + e2*sin(phi) - (b/2)*sin(gamma-beta);
C_y2 = A_y + e2*cos(phi) + (b/2)*cos(gamma-beta);

%% Plots

% Camber angle vs Deflection angle
theta2_deg = theta2*(180/pi);
camber_deg = camber*(180/pi); 
figure(1)
plot(theta2_deg, camber_deg)
hold off
ylabel("Camber Angle, deg");
xlabel("Lower Wishbone Angle, deg");

if answer2 == "real-world"
    % Camber angle vs Chassis Vertical Displacement
    camber_deg = camber*(180/pi);
    figure(2)
    plot(delta_x*1000, camber_deg)
    hold off
    ylabel("Camber Angle, deg");
    xlabel("Chassis Vertical Displacement, mm");
end

% Stub axle vs Deflection angle
figure(3)
yyaxis left
plot(theta2_deg, C_x)
hold off
ylabel("Stub axle x position");
yyaxis right
plot(theta2_deg, C_y)
ylabel("Stub axle y position");
xlabel("Lower Wishbone Angle");

% Animation
figure(4)
for i = 1:length(theta2)
    plot(A_x(i),A_y(i),'o', M_x(i),M_y(i),'o', B_x(i),B_y(i),'o', N_x(i),N_y(i),'o', C_x(i),C_y(i),'o', C_x2(i),C_y2(i),'o', 'Color', '#035EAB')
    hold on

    % Line AM
    x1 = [A_x(i), M_x(i)];
    y1 = [A_y(i), M_y(i)];
    plot(x1,y1,'Color', '#035EAB')
    hold on
    
    % Line BN
    x2 = [B_x(i), N_x(i)];
    y2 = [B_y(i), N_y(i)];
    plot(x2,y2, 'Color', '#035EAB')
    hold on

    % Line C1C2
    x3 = [C_x(i), C_x2(i)];
    y3 = [C_y(i), C_y2(i)];
    plot(x3,y3, 'Color', '#035EAB')
    hold on

    % Line AB
    x4 = [A_x(i), B_x(i)];
    y4 = [A_y(i), B_y(i)];
    plot(x4,y4, 'Color', '#035EAB')
    hold on

    % Coordinates of wheel corners on animation
    % top left corner
    tl_x(i) = C_x2(i) - tyre_width/4;
    tl_y(i) = C_y2(i) + tyre_height/2;

    % top right corner
    tr_x(i) = C_x2(i) + 3*tyre_width/4;
    tr_y(i) = C_y2(i) + tyre_height/2;

    % bottom left corner
    bl_x(i) = C_x2(i) - tyre_width/4;
    bl_y(i) = C_y2(i) - tyre_height/2;

    % bottom right corner
    br_x(i) = C_x2(i) + 3*tyre_width/4;
    br_y(i) = C_y2(i) - tyre_height/2;

    hold on

    % Plot the wheel
    pgon = polyshape([bl_x(i) br_x(i) tr_x(i) tl_x(i)],[bl_y(i) br_y(i) tr_y(i) tl_y(i)]);
    plot(rotate(pgon, -camber_deg(i), [C_x2(i) C_y2(i)]), 'FaceColor', '#035EAB');
    hold off

    axis equal
    ylabel("Y Position, mm");
    xlabel("X Position, mm");
    %pause(0.1) % add a pause if you'd like to watch the animation more slowly
    drawnow
end