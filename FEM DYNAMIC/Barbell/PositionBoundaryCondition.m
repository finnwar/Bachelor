% Implementation of Boundary Conditions

function [U_Boundary, affectedNodes, U_dot,U_ddot] = PositionBoundaryCondition(NodeGrid, t)
    %Initialization  
    U_Boundary = zeros(NodeGrid(end, end),1);
    U_dot =zeros(size(U_Boundary));
    U_ddot = zeros(size(U_Boundary));
    % Motion function of forced movement
    
    amplitude = 0.1; %[m]
    frequency = 1;   %[Hz]
    phase = 0;  %[s]

    Position = @(t) amplitude*sin(2*pi*frequency*t-phase);
    Velocity = @(t) amplitude*2*pi*frequency*cos(2*pi*frequency*t-phase);
    Acceleration = @(t) -amplitude*4*pi^2*frequency^2*sin(2*pi*frequency*t-phase);

    % Left side x-Displacement
    
    U_Boundary(NodeGrid(1,:)) = Position(t);
    U_dot(NodeGrid(1,:)) = Velocity(t);
    U_ddot(NodeGrid(1,:)) = Acceleration(t);
    % Left side y displacement
    
    U_Boundary(NodeGrid(2,:)) = 0;
    U_dot(NodeGrid(2,:)) = 0;
    U_ddot(NodeGrid(2,:)) = 0;
    
    
    affectedNodes = [NodeGrid(1,:) NodeGrid(2,:)];
    

end