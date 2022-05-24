clear all
close all
clc

%   ______________________________________________________________________
%   **********************************************************************
%   Script uses elements proposed in:
%   MMehrez,MPC-and-MHE-implementation-in-MATLAB-using-Casadi,(2021),GitHub
%   https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%
%   Description:
%   This script is used to demonstrate the  stablity of the controller.
%   subject to impulse input.
%   ______________________________________________________________________
%   **********************************************************************


% CasADi v3.5.5
% To use code the addpath will have to change to the same location the
% folder was downloaded under
addpath('C:\Users\phili\OneDrive - University of Bath\Year 3\Semester 2\Matlab\MPC_Controller\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

h = 0.05;   % time step [s]
N = 40;     % prediction horizon
car_diam = 1; % diamter of car [m]
laps = 1;   % number of laps

% cone coordinates x dim
conesx= [1; 1; 2; 2; 3; 3; 4; 4; 5; 5; 6; 6; 7; 7; 8; 8; 9; 9; 10; 10;...
    11; 11; 12; 12; 13; 13; 14; 14; 15; 15; 16; 16; 17; 17; 18; 18; 19; 19;...
    20; 20; 21; 21; 22; 22; 23; 23; 24; 24; 25; 25; 26; 26; 27; 27; 28; 28; 29; 29];
% cone coordinates y dim
conesy= [10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10;...
    -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10;...
    10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10; 10;...
    -10; 10; -10; 10; -10; 10; -10; 10; -10; 10; -10];

% reformulate for number of laps
conesx = repmat(conesx, laps, 1);  
conesy = repmat(conesy, laps, 1);  

cones_diam = 0.5; % meters
cone_horizon = 6; % number of cones

x_ref = [[1 ; 0; 0; 1],...
         [2 ; 0; 0; 1],...
         [3 ; 0; 0; 1],...
         [4 ; 1000; 0; 1000],...
         [5 ; 0; 0; 1],...
         [6 ; 0; 0; 1],...
         [7 ; 0; 0; 1],...
         [8 ; 0; 0; 1],...
         [9 ; 0; 0; 1],...
         [10; 0; 0; 1],...
         [11; 0; 0; 1],...
         [12; 0; 0; 1],...
         [13; 0; 0; 1],...
         [14; 0; 0; 1],...
         [15; 0; 0; 1],...
         [16; 0; 0; 1],...          
         [17; 0; 0; 1],...
         [18; 0; 0; 1],...
         [19; 0; 0; 1],...
         [20; 0; 0; 1],...
         [21; 0; 0; 1],...
         [22; 0; 0; 1],...
         [23; 0; 0; 1],...
         [24; 0; 0; 1],...
         [25; 0; 0; 1],...
         [26; 0; 0; 1],...
         [27; 0; 0; 1],...
         [28; 0; 0; 1],...
         [29; 0; 0; 1]];
     
          % reference car positions and velocities
        
x_ref = repmat(x_ref, 1, laps); %reformaulate for number of laps

% retrieve and store car paramenters
parameters = Variables();
m = parameters.m;
Lr = parameters.Lr;
Lf = parameters.Lf;
Iz = parameters.Iz;
Fz = parameters.Frz*1000;

% Control input max/min values
Frx_max = 500; Frx_min = -500;
sigma_max = 0.55; sigma_min = -sigma_max;

% define symbolic states
x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
vx = SX.sym('vx'); vy = SX.sym('vy');  omega = SX.sym('omega');
states = [x;y;theta;vx;vy;omega]; n_states = length(states);

% define symbolic control inputs
Frx = SX.sym('Frx'); sigma = SX.sym('sigma');
controls = [Frx;sigma]; n_controls = length(controls);

% Slip Angle Calculations
alphaf = -atan((omega*Lf + vy)/vx)+(sigma);
alphar = atan((omega*Lr - vy)/vx);
 
% Laterial Tyre Force Calculation
Ffy = 200*alphaf;
Fry = 200*alphar;

% non-linear kinematic bicycle model
rhs = [vx*cos(theta)-vy*sin(theta);...
    vx*sin(theta)+vy*cos(theta);...
    omega;...
    (Frx-Ffy*sin(sigma))/m+vy*omega;...
    (1/m)*(Fry+Ffy*cos(sigma)-vx*omega);...
    (1/Iz)*(Ffy*cos(sigma)*Lf-Fry*Lr)]; 

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)

% parameters (which include the initial state and the reference state)
P = SX.sym('P',n_states + n_states-2);

% A vector that represents the states over the optimization problem.
X = SX.sym('X',n_states,(N+1));

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(4,4); Q(1,1) = 20; Q(2,2) = 20; Q(3,3)= 0.01;...
    Q(4,4) = -10; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.001; R(2,2) = 40; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:6)]; % initial condition constraints

% RK-4 with Multiple Shooting method 
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st(1:4)-P(7:10))'*Q*(st(1:4)-P(7:10)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    k1 = f(st, con);   
    k2 = f(st + h/2*k1, con); 
    k3 = f(st + h/2*k2, con); 
    k4 = f(st + h*k3, con); 
    st_next_RK4=st +h/6*(k1 +2*k2 +2*k3 +k4);    
    g = [g;st_next-st_next_RK4]; % compute constraints 
end

OPT_variables = [reshape(X,6*(N+1),1);reshape(U,2*N,1)];

% Casadi Solver Attributes
opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;


args = struct;

args.lbg(1:6*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:6*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbg(6*(N+1)+1 : 6*(N+1)+cone_horizon*(N+1)) = -inf; %Track constraints
args.ubg(6*(N+1)+1 : 6*(N+1)+cone_horizon*(N+1)) = 0; %Track constraints  

args.lbx(1:6:6*(N+1),1) = -200; %state x lower bound
args.ubx(1:6:6*(N+1),1) = 200; %state x upper bound
args.lbx(2:6:6*(N+1),1) = -100; %state y lower bound
args.ubx(2:6:6*(N+1),1) = 100; %state y upper bound
args.lbx(3:6:6*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:6:6*(N+1),1) = inf; %state theta upper bound
args.lbx(4:6:6*(N+1),1) = 1; %state vx lower bound
args.ubx(4:6:6*(N+1),1) = 5; %state vx upper bound
args.lbx(5:6:6*(N+1),1) = -inf; %state vy lower bound
args.ubx(5:6:6*(N+1),1) = inf; %state vy upper bound
args.lbx(6:6:6*(N+1),1) = -inf; %state omega lower bound
args.ubx(6:6:6*(N+1),1) = inf; %state omega upper bound

args.lbx(6*(N+1)+1:2:6*(N+1)+2*N,1) = Frx_min; %Frx lower bound
args.ubx(6*(N+1)+1:2:6*(N+1)+2*N,1) = Frx_max; %Frx upper bound
args.lbx(6*(N+1)+2:2:6*(N+1)+2*N,1) = sigma_min; %sigma lower bound
args.ubx(6*(N+1)+2:2:6*(N+1)+2*N,1) = sigma_max; %sigma upper bound

%----------------------------------------------
% ABOVE IS PROBLEM SET UP

t0 = 0;
x0 = [0 ; 0 ; 0 ;  1; 0 ; 0];    % initial condition.
u0 = zeros(N,2);        % initialization of the control inputs
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables
 
sim_tim = 40; % Maximum simulation time

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;
g0 = g;

mpc_iter = 0; 
xx1 = [];
u_cl=[];

reference_pos = cone_horizon/2; % counter for car reference position
cone_sample = 1;                % counter for car 

main_loop = tic;

while(mpc_iter < 40 / h && cone_sample+cone_horizon-1<length(conesx))
    
    cones_x = conesx(cone_sample:cone_sample+cone_horizon-1);
    cones_y = conesy(cone_sample:cone_sample+cone_horizon-1);   

    xs = x_ref(:,reference_pos);

    % Add constraints for cone avoidance
    g=g0;
    for k = 1:N+1   
        g = [g ; -sqrt((X(1,k)-cones_x).^2+(X(2,k)-cones_y).^2) + (car_diam/2 + cones_diam/2)];
    end
    
    % set optimization solver parameters
    nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
    solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

    % finds line which car needs to cross to update cones and reference
    grad = (cones_y(3)-cones_y(4))/(cones_x(3)-cones_x(4));
    if (abs(grad) == inf)
        c = cones_x(3);
    elseif (grad == 0)
        c = cones_y(3); 
    else
        c = cones_y(3)-grad*cones_x(3); 
    end
    
    trigger = 0;
    % Start MPC controller to optimize control inputs to given reference
    % exits if the controller takes too long or if the controller needs
    % updating
    while(trigger==0 && mpc_iter < 40 / h)
        
        % exits mpc controller if car has passed second set of cones in
        % sample
        if grad == 0
            if (cones_x(4)>cones_x(3) && x0(2)<c)
                trigger = 1;
            elseif (cones_x(4)<cones_x(3) && x0(2)>c)
                trigger = 1;
            end
        elseif (grad == inf && x0(1)>c)
            trigger = 1;
        elseif (grad == -inf && x0(1)<c)
            trigger = 1;
        elseif (abs(grad) ~= inf && grad ~= 0)
            if(cones_y(3)>cones_y(4))
                if (grad>0 && x0(2)<grad*x0(1)+c)
                    trigger = 1;
                elseif (grad<0 && x0(2)>grad*x0(1)+c)
                    trigger = 1;
                end
            else
                if (grad>0 && x0(2)>grad*x0(1)+c)
                    trigger = 1;
                elseif (grad<0 && x0(2)<grad*x0(1)+c)
                    trigger = 1;
                end
            end
        end

        args.p  = [x0;xs]; % set the values of the parameters vector
        args.x0  = [reshape(X0',6*(N+1),1);reshape(u0',2*N,1)]; % initial value of the optimization variables
         sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        u = reshape(full(sol.x(6*(N+1)+1:end))',2,N)';  % get controls only from the solution
        xx1(:,1:6,mpc_iter+1)= reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get prediction 
        u_cl= [u_cl ; u(1,:)];
        t(mpc_iter+1) = t0;
        [t0, x0, u0] = shift(h, t0, x0, u,f); % Apply the control and shift the solution
        xx(:,mpc_iter+2) = x0;
        X0 = reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution 
        X0 = [X0(2:end,:);X0(end,:)];        % Shift trajectory to initialize the next step
        mpc_iter = mpc_iter + 1;

    end
    reference_pos=reference_pos+1;
    cone_sample = cone_sample+2;

end

main_loop_time = toc(main_loop);
average_mpc_time = main_loop_time/(mpc_iter+1)

Draw_MPC_impulse(t,xx,xx1,u_cl,x_ref,N,car_diam,conesx,conesy,cones_diam,cone_horizon)