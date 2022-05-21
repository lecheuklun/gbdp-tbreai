 %% Load the data
% First load the data into the workspace
% For this example, the data is stored as a Table and has been already
% filtered, cropped and pre-processed.
% The data in this example is already in ISO-W and all channels in SI units
% (N, Nm, m, s, rad, Pa)

clear; clc;
% load('TyreData.mat')
load('TyreDataRear.mat')

plot(TyreData.P)

%filter data
TyreData = TyreData(1:10:14000,:);

%% Fitting process

% Load a TIR file as a starting point for the fitting
InitialParameterSet = mfeval.readTIR('PacejkaBook_Defaults.tir');

% Set nominal parameters of the model (DO NOT CHANGE AFTER)
InitialParameterSet.UNLOADED_RADIUS = 0.26; % Unloaded tire radius
InitialParameterSet.FNOMIN = 1200; % Nominal load
InitialParameterSet.LONGVL = 16.7; % Nominal reference speed
InitialParameterSet.NOMPRES = 85000; % Nominal inflation pressure

% Create the initial parameters for the fitting (seeds)
x0(1) = InitialParameterSet.PCY1;
x0(2) = InitialParameterSet.PDY1;
x0(3) = InitialParameterSet.PDY2;
x0(4) = InitialParameterSet.PDY3;
x0(5) = InitialParameterSet.PEY1;
x0(6) = InitialParameterSet.PEY2;
x0(7) = InitialParameterSet.PEY3;
x0(8) = InitialParameterSet.PEY4;
x0(9) = InitialParameterSet.PEY5;
x0(10) = InitialParameterSet.PKY1;
x0(11) = InitialParameterSet.PKY2;
x0(12) = InitialParameterSet.PKY3;
x0(13) = InitialParameterSet.PKY4;
x0(14) = InitialParameterSet.PKY5;
x0(15) = InitialParameterSet.PKY6;
x0(16) = InitialParameterSet.PKY7;
x0(17) = InitialParameterSet.PHY1;
x0(18) = InitialParameterSet.PHY2;
x0(19) = InitialParameterSet.PVY1;
x0(20) = InitialParameterSet.PVY2;
x0(21) = InitialParameterSet.PVY3;
x0(22) = InitialParameterSet.PVY4;
x0(23) = InitialParameterSet.PPY1;   %influence of inflation pressure on cornering stiffness
x0(24) = InitialParameterSet.PPY2;   %influence of inflation pressure on dependency of nominal tyre load on cornering stiffness
x0(25) = InitialParameterSet.PPY3;   %linear influence of inflation pressure on lateral peak friction
x0(26) = InitialParameterSet.PPY4;   %quadratic influence of inflation pressure on lateral peak friction
x0(27) = InitialParameterSet.PPY5;

% Declare the anonymous function (Cost function) for the fitting
% The @ operator creates the handle, and the parentheses () immediately
% after the @ operator include the function input arguments
fun = @(X) costFyPure(X, TyreData, InitialParameterSet);

% Options for the fitting function lsqnonlin
options.TolFun = 1e-08; % Low tolerance to ensure good fitting
options.MaxFunEvals = 999999; % Very high to avoid this stop criteria
options.MaxIter = 99999;
options.Display = 'iter';

% Non-linear least squares fitting formula
% lsqnonlin will try to minimize the output of the cost function (error).
% Go to the cost function "costFyPure" to check how this is performed
tic
X_OPTIM = lsqnonlin(fun,x0,[],[],options);
toc

% Create a copy of the initial parameters and replace the fitted parameters
OptimParameterSet = InitialParameterSet;

OptimParameterSet.PCY1 = X_OPTIM(1);
OptimParameterSet.PDY1 = X_OPTIM(2);
OptimParameterSet.PDY2 = X_OPTIM(3);
OptimParameterSet.PDY3 = X_OPTIM(4);
OptimParameterSet.PEY1 = X_OPTIM(5);
OptimParameterSet.PEY2 = X_OPTIM(6);
OptimParameterSet.PEY3 = X_OPTIM(7);
OptimParameterSet.PEY4 = X_OPTIM(8);
OptimParameterSet.PEY5 = X_OPTIM(9);
OptimParameterSet.PKY1 = X_OPTIM(10);
OptimParameterSet.PKY2 = X_OPTIM(11);
OptimParameterSet.PKY3 = X_OPTIM(12);
OptimParameterSet.PKY4 = X_OPTIM(13);
OptimParameterSet.PKY5 = X_OPTIM(14);
OptimParameterSet.PKY6 = X_OPTIM(15);
OptimParameterSet.PKY7 = X_OPTIM(16);
OptimParameterSet.PHY1 = X_OPTIM(17);
OptimParameterSet.PHY2 = X_OPTIM(18);
OptimParameterSet.PVY1 = X_OPTIM(19);
OptimParameterSet.PVY2 = X_OPTIM(20);
OptimParameterSet.PVY3 = X_OPTIM(21);
OptimParameterSet.PVY4 = X_OPTIM(22);
OptimParameterSet.PPY1 = X_OPTIM(23);
OptimParameterSet.PPY2 = X_OPTIM(24);
OptimParameterSet.PPY3 = X_OPTIM(25);
OptimParameterSet.PPY4 = X_OPTIM(26);
OptimParameterSet.PPY5 = X_OPTIM(27);

%% Plot results

% Filter data to plot specific conditions:
indFz1 = TyreData.Fz > 1000 & TyreData.Fz < 1200;   % 1100 N
indFz2 = TyreData.Fz > 825 & TyreData.Fz < 925;   % 875 N
indFz3 = TyreData.Fz > 600 & TyreData.Fz < 700;    % 650 N
indFz4 = TyreData.Fz > 400 & TyreData.Fz < 500;    % 450 N
indFz5 = TyreData.Fz > 175 & TyreData.Fz < 275;    % 225 N
indIA = TyreData.IA > -0.01 & TyreData.IA < 0.01;   % 0 rad
indP = TyreData.P > 8e4 & TyreData.P < 9e4;
% indP = TyreData.P > 80 & TyreData.P < 90;% 83160 Pa
indFz = indFz1 | indFz2 | indFz3 | indFz4 | indFz5;
filt = indFz & indIA & indP;

% Create data inputs to do a data replay with MFeval and check the fitting
% quality
evalFz1 = ones(100,1)*1100;
evalFz2 = ones(100,1)*875;
evalFz3 = ones(100,1)*650;
evalFz4 = ones(100,1)*450;
evalFz5 = ones(100,1)*225;
evalNull = zeros(100, 1);
evalSA = linspace(-1,1)'; %changed
evalVx = ones(100, 1)*16;
evalP = ones(100,1)*83160;

MFinput1 = [evalFz1 evalNull evalSA evalNull evalNull evalVx evalP];
MFinput2 = [evalFz2 evalNull evalSA evalNull evalNull evalVx evalP];
MFinput3 = [evalFz3 evalNull evalSA evalNull evalNull evalVx evalP];
MFinput4 = [evalFz4 evalNull evalSA evalNull evalNull evalVx evalP];
MFinput5 = [evalFz5 evalNull evalSA evalNull evalNull evalVx evalP];

% Call mfeval with the optimized parameters
MFout1 = mfeval(OptimParameterSet,MFinput1,111);
MFout2 = mfeval(OptimParameterSet,MFinput2,111);
MFout3 = mfeval(OptimParameterSet,MFinput3,111);
MFout4 = mfeval(OptimParameterSet,MFinput4,111);
MFout5 = mfeval(OptimParameterSet,MFinput5,111);

% Plot data vs Fitted Model
figure
hold on
plot(TyreData.SA(filt), TyreData.Fy(filt),'o')
plot(MFout1(:,8), MFout1(:,2),'-', 'linewidth', 2)
plot(MFout2(:,8), MFout2(:,2),'-', 'linewidth', 2)
plot(MFout3(:,8), MFout3(:,2),'-', 'linewidth', 2)
plot(MFout4(:,8), MFout4(:,2),'-', 'linewidth', 2)
plot(MFout5(:,8), MFout5(:,2),'-', 'linewidth', 2)
grid on
xlabel('Slip Angle [rad]')
ylabel('Lateral Force [N]')
title('PureFy fitting')
legend('Data', 'Model: Fz=1200N', 'Model: Fz=1650N', 'Model: Fz= 760N', 'Model: Fz= 450N', 'Model: Fz= 225N')


%% Nested functions
function [ error ] = costFyPure(X, Data, ParameterSet)
%COSTFYPURE calls MFeval and calculates the error between the model and the
%input data.
%
% error = costFyPure(X, Data, ParameterSet)
%
% X: Is a structure that contains the FyPure parameters that are being
%       fitted. X is changing all the time when lsqnonlin is calling this
%       function.
% Data: Is a Table that contains the Data used to measure the error
%       of the model that is being fitted.
% ParameterSet: Is a structure of MF6.1 parameters. The parameters are used
%       only to call MFeval without errors.
%
% Example:
% error = costFyPure(Xstructure, TableData, ParameterSet)

% Create the Inputs for MFeval
INPUTS = [Data.Fz Data.SR Data.SA Data.IA Data.Phit Data.Vx Data.P Data.W];

% Select use mode 211. For more info go to the documentation of MFeval
USE_MODE = 211;

% Unpack the parameters that are being fitted and replace them into the
% ParameterSet.
ParameterSet.PCY1 = X(1);
ParameterSet.PDY1 = X(2);
ParameterSet.PDY2 = X(3);
ParameterSet.PDY3 = X(4);
ParameterSet.PEY1 = X(5);
ParameterSet.PEY2 = X(6);
ParameterSet.PEY3 = X(7);
ParameterSet.PEY4 = X(8);
ParameterSet.PEY5 = X(9);
ParameterSet.PKY1 = X(10);
ParameterSet.PKY2 = X(11);
ParameterSet.PKY3 = X(12);
ParameterSet.PKY4 = X(13);
ParameterSet.PKY5 = X(14);
ParameterSet.PKY6 = X(15);
ParameterSet.PKY7 = X(16);
ParameterSet.PHY1 = X(17);
ParameterSet.PHY2 = X(18);
ParameterSet.PVY1 = X(19);
ParameterSet.PVY2 = X(20);
ParameterSet.PVY3 = X(21);
ParameterSet.PVY4 = X(22);
ParameterSet.PPY1	=  X(23)   	;%influence of inflation pressure on cornering stiffness
ParameterSet.PPY2	=  X(24)  	;%influence of inflation pressure on dependency of nominal tyre load on cornering stiffness
ParameterSet.PPY3	=  X(25)   	;%linear influence of inflation pressure on lateral peak friction
ParameterSet.PPY4	=  X(26)   	;%quadratic influence of inflation pressure on lateral peak friction
ParameterSet.PPY5	=  X(27)   	;

% Call MFeval
OUTPUT = mfeval(ParameterSet,INPUTS,USE_MODE);

% Get the Fy from the MF6.1 model
Fy_MFeval = OUTPUT(:,2);

% Calculate error against the data
error = (Data.Fy - Fy_MFeval);
end
