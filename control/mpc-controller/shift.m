%   ______________________________________________________________________
%   **********************************************************************
%   Taken from:
%   MMehrez,MPC-and-MHE-implementation-in-MATLAB-using-Casadi,(2021),GitHub
%   https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%   ______________________________________________________________________
%   **********************************************************************


function [t0, x0, u0] = shift(T, t0, x0, u,f)
%This function takes the calculated control action and current state of the
%vehicle and applies the control action using an euler approximation. The
%new location is then returned and the time step is updated along with the
%updated predicted control action.

st = x0;
con = u(1,:)';
f_value = f(st,con);
st = st+ (T*f_value);
x0 = full(st);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end
