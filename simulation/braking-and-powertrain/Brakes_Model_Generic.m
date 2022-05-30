function Brakes_Model_Generic(t,dt,Front_Valve_Position,Rear_Valve_Position,P_High,P_Calip,P_Low)

X(1) = Front_Valve_Position;
X(2) = Rear_Valve_Position;

for i = 1:2

    if X(i) > 0
        P_Ref = P_High;
    else
        P_Ref = P_Low;
    end
    

    if abs(P_Ref - P_Calip) > 0.1
        
        dPA = State_Deriv_Calip_P(t,P_Calip,X(i),P_Ref);
        A = dPA*dt;
        
        dPB = State_Deriv_Calip_P(t+dt/2,P_Calip+A/2,X(i),P_Ref);
        B = dPB*dt;
        
        dPC = State_Deriv_Calip_P(t+dt/2,P_Calip+B/2,X(i),P_Ref);
        C = dPC*dt;
        
        dPD = State_Deriv_Calip_P(t+dt,P_Calip+C,X(i),P_Ref);
        D = dPD*dt;
        
        Braking_Pressure(i) = P_Calip + (1/6)*(A + 2*B + 2*C + D)
        
    end
end
    Caliper_Piston_Area = pi*((25.4*10^-3)/2)^2;
    Pad_Coefficient_of_Friction = 0.39; %APH420 coefficient of friction
    %Pad_Coefficient_of_Friction = 0.41; %RQ3 coefficient of friction
    
   Front_Braking_Torque = Braking_Pressure(1)*10^5*Caliper_Piston_Area*8*Pad_Coefficient_of_Friction*0.08;
   Rear_Braking_Torque = Braking_Pressure(2)*10^5*Caliper_Piston_Area*4*Pad_Coefficient_of_Friction*0.08;
    
   Total_Longditudinal_Braking_Force = (Front_Braking_Torque + Rear_Braking_Torque)/0.2
    
   t = t+dt;