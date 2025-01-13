function [f_cond,affectedNodes] = ForceBoundaryCondition(NodeGrid,t)
    
    %Initialization

    f_cond = zeros(NodeGrid(end,end),1);
    
    % Condition

    maxForce = 0.0001; %[N(/m)]
    frequency = 10; %[Hz]
    phase = 0; %[-]

    forceFunction = @(t) maxForce * cos(2*pi*frequency*t - phase);

    % Application of force at the top middle of part (singe node!!!)

    
    node = NodeGrid(end/2);

    if mode(node,2)==1
        f_cond(node) = forceFunction(t);
        affectedNodes = node;
    else
        f_cond(node+1) = forceFunction(t);
        affectedNodes = node+1;
    end



end