function [f_cond,affectedNodes] = ForceBoundaryCondition(NodeGrid,t)
    
    %Initialization

    f_cond = zeros(NodeGrid(end,end),1);
    
    % Condition

    maxForce = 1000; %[N]
    frequency = 10; %[Hz]
    phase = 0; %[-]

    forceFunction = maxForce * cos(2*pi*frequency*t - phase);

    % Application of force at the top middle of part (singe node!!!)

    
    node = NodeGrid(end/2, end);

    if mod(node,2)==1
        f_cond(node) = forceFunction;
        affectedNodes = node;
    else
        f_cond(node-1) = forceFunction;
        affectedNodes = node-1;
    end



end