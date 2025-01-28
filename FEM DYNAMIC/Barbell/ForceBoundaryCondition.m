function [f_cond,node] = ForceBoundaryCondition(NodeGrid,t)
    
    %Initialization

    f_cond = zeros(NodeGrid(end,end),1);
    
    % Condition

    maxForce = -100; %[N]
    frequency = 1/8; %[Hz]
    phase = 0; %[-]

    forceFunction = maxForce * sin(2*pi*frequency*t - phase);

    % Application of force at the top end of part (single node!!!)

    
    node = NodeGrid(end-1, end);

    if mod(node,2)==1
        f_cond(node) = forceFunction;
    else
        f_cond(node-1) = forceFunction;
        node = node-1;
    end



end