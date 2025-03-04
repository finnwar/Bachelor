function [f_cond,node] = ForceBoundaryCondition(NodeGrid,t)
    
    %Initialization

    f_cond = zeros(NodeGrid(end,end),1);
    
    % Condition

    maxForce = 100; %[N]
    frequency = 1/8; %[Hz]
    phase = 0; %[-]

    forceFunction = maxForce * sin(2*pi*frequency*t - phase);

    % Application of force at the top end of part (single node!!!)

    
    node = NodeGrid(end, end);

    f_cond(node) = forceFunction;
end