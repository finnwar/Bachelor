function [U] = StaticFEM (K,M,NodeGrid)
    
    [U_Boundary, ConstrainedNodes] = PositionBoundaryCondition(NodeGrid,0);
    [f,~] = ForceBoundaryCondition(NodeGrid,0);
    
    f_tilde = f-K*U_Boundary;
    
    f_tilde(ConstrainedNodes) = [];
    
    %Reduce system
    K_tilde = K;
    K_tilde(:,ConstrainedNodes) = [];
    K_tilde(ConstrainedNodes,:) = [];

    

    U_tilde = K_tilde\f_tilde;
    U=U_tilde;
    
    %Reinsert the conditional displacements into displacement vector.
    ConstrainedNodes=sort(ConstrainedNodes);  %The vector containing the boundary nodes has to be sorted for the following method

    for i = 1:length(ConstrainedNodes)
       U = [U(1:(ConstrainedNodes(i)-1)); U_Boundary(ConstrainedNodes(i)); U((ConstrainedNodes(i)):length(U))];
    end

end
