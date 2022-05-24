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
%   This script is used to demonstrate the  capabilities of the controller
%   around a simulated track with simulated reference inputs. This
%   controller uses a RK-4 approximation of the Dynamic Bicycle model of
%   the vehicle with linear tire model.
%   ______________________________________________________________________
%   **********************************************************************


% CasADi v3.5.5
% To use code the addpath will have to change to the same location the
% folder was downloaded under
addpath('C:\Users\phili\OneDrive - University of Bath\Year 3\Semester 2\Matlab\MPC_Controller\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

h = 0.05;   % time step [s]
N = 60;     % prediction horizon
car_diam = 1; % diamter of car [m]
laps = 1;   % number of laps

% cone coordinates x dim
conesx=[2;2;5;5;8;8;9.5;12;10.5;14;12;15;14.5;16.5;17;18;19;19.5;21;22.5;...
    21.5;24;21.5;25;20;23.5;19.5;20.5;18;17;15;15;12;12;...
    10;10;8;8;6;6;4.5;4.5;3;3;2;1;1;-2;0;-4;0;-4.5;1;-1]; 

% cone coordinates y dim
conesy=[1;-2;1;-2;1;-2;1.5;-1;2.5;0.5;3.5;2;4.5;2;5;2;4.5;2;5;2;...
    6;4;7;6.5;7.5;8.5;7.5;10;7;9.5;6.5;9;6;8.5;...
    6;8;6;8;6;8;5;8;5;8;5;8;4;7.5;3;4.5;1.5;0.5;1;-2];    

% reformulate for number of laps
conesx = repmat(conesx, laps, 1);  
conesy = repmat(conesy, laps, 1);  

cones_diam = 0.5; % meters
cone_horizon = 8; % number of cones
v_ref = 4;        % reference velocity

x_ref = [2, -0.5;
        5 , -0.5;
        8 , -0.5;
        11 , 0.5; 
        12 , 1.5;
        13.5 , 2.75; 
        15 , 3.25; 
        17 , 3.5; 
        19 , 3; 
        22 , 3.5; 
        23 , 5; 
        23 , 6.75; 
        22 , 8; 
        20 , 9;
        17.5 , 8.5;
        15 , 7.75; 
        12 , 7.25; 
        10 , 7; 
        8 , 7; 
        6 , 7;
        4.5 , 6.5; 
        3 , 6.5;
        1.5, 6  ;
        -1 , 5.5; 
        -2 , 3.5; 
        -2 , 1; 
        0 , -0.5];;  % reference car positions and velocities
    
   
x_ref = positions_to_ref(x_ref, v_ref); %conversiion of waypoints to refernec states  
x_ref = repmat(x_ref, 1, laps); %reformaulate for number of laps

% retrieve and store car paramenters
parameters = Variables();
m = parameters.m;
Lr = parameters.Lr;
Lf = parameters.Lf;
Iz = parameters.Iz;
Fz = parameters.Frz*1000;

% Control input max/min values
Frx_max = 3000; Frx_min = -3630;
sigma_max = 0.55; sigma_min = -sigma_max;

% define symbolic states
x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
vx = SX.sym('vx'); vy = SX.sym('vy');  omega = SX.sym('omega');
states = [x;y;theta;vx;vy;omega]; n_states = length(states);

% define symbolic control inputs
Frx = SX.sym('Frx'); sigma = SX.sym('sigma');
controls = [Frx;sigma]; n_controls = length(controls);

% Slip Angle Calculations
alphaf = -atan((omega*Lf + vy)/(vx+0.1))+(sigma);
alphar = atan((omega*Lr - vy)/(vx+0.1));
 
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
    Q(4,4) = 1; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.0001; R(2,2) = 55; % weighing matrices (controls)

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
args.ubx(4:6:6*(N+1),1) = 10; %state vx upper bound
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
x0 = [0 ; 0 ; 0 ; 0; 0 ; 0];    % initial condition.
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
cone_sample = 1;                % counter for cone positions

main_loop = tic;

while(reference_pos<=length(x_ref))
    
    %Sample horizon of cones to predict over    
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
    solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

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
    while(trigger==0 && mpc_iter < 40*laps / h)
        
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
        xx(:,mpc_iter+2) = x0;                % Store location of vehicle
        X0 = reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution 
        X0 = [X0(2:end,:);X0(end,:)];        % Shift trajectory to initialize the next step
        mpc_iter = mpc_iter + 1;

    end
    reference_pos = reference_pos+1;
    cone_sample = cone_sample+2;

end

main_loop_time = toc(main_loop);
average_mpc_time = main_loop_time/(mpc_iter+1)

Draw_MPC_track(t,xx,xx1,u_cl,x_ref,N,car_diam,conesx,conesy,cones_diam,cone_horizon)