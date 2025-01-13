% Implementation of Boundary Conditions

function [U_Boundary, affectedNodes] = PositionBoundaryCondition(NodeGrid, t)
    %Initialization  
    U_Boundary = zeros(NodeGrid(end, end),1);
    
    % Motion function of forced movement
    
    amplitude = 0.1; %[m]
    frequency = 1;   %[Hz]
    phase = 0;  %[s]

    PosFunction = @(t) amplitude*sin(2*pi*frequency*t-phase);


    % Left side x-Displacement

    U_Boundary(NodeGrid(1,:)) = PosFunction(t);

    % Left side y displacement
    
    U_Boundary(NodeGrid(2,:)) = 0;
    
    affectedNodes = [NodeGrid(1,:) NodeGrid(2,:)];
    

end