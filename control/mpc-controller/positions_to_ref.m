function [Pos] = positions_to_ref(Positions, Set_velocity)
%PositionToReference converts x-y coordinates of path to reference state
%   X-Y coordinates are used to calculate the yaw angle required for the
%   car to follow the track. The positions have been calcualted by the path
%   planning team. This data is then combined with a value for the
%   reference velocity for a full set of states for the MPC controller

v_ref = ones(length(Positions), 1).*Set_velocity;

for n = 1:length(Positions(:,1))-1
    
    
    Point1 = Positions(n,:);
    Point2 = Positions(n+1,:);
    
    % Determine spline from each adjacent set of track waypoints
    X = Point2(1)-Point1(1);
    Y = Point2(2)-Point1(2);
    
    % Calculate angle of splie from coordinate axis
    YawAngle(n) = atan2(Y,X);

end

Pos = [Positions, [YawAngle'; 0], v_ref]';

