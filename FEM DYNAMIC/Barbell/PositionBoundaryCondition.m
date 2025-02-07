% Implementation of Boundary Conditions

function [U_Boundary, affectedNodes, U_dot,U_ddot] = PositionBoundaryCondition(NodeGrid, t)
    %Initialization  
    U_Boundary = zeros(NodeGrid(end, end),1);
    U_dot =zeros(size(U_Boundary));
    U_ddot = zeros(size(U_Boundary));
    % Motion function of forced movement
    
    amplitude = 0.01; %[m]
    frequency = 0.001;   %[Hz]
    phase = 0;  %[-]

    Position = amplitude*sin(2*pi*frequency*t-phase);
    Velocity = amplitude*2*pi*frequency*cos(2*pi*frequency*t-phase);
    Acceleration =amplitude*4*pi^2*frequency^2*sin(2*pi*frequency*t-phase);


    for i = 1:length(NodeGrid(1,:))
        % Left side x-Displacement
        U_Boundary(NodeGrid(1,i)) = Position;
        U_dot(NodeGrid(1,i)) = Velocity;
        U_ddot(NodeGrid(1,i)) = Acceleration;
        
        % Left side y displacement
        U_Boundary(NodeGrid(2,i)) = 0;
        U_dot(NodeGrid(2,i)) = 0;
        U_ddot(NodeGrid(2,i)) = 0;
    end
    
    affectedNodes = [NodeGrid(1,:) NodeGrid(2,:)];
    

end