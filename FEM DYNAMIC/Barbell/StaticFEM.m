function [U] = StaticFEM (K,M,NodeGrid)
    
    [U_Boundary, BoundaryNodes] = PositionBoundaryCondition(NodeGrid,0);
    [f,~] = ForceBoundaryCondition(NodeGrid,0);
    
    f_tilde = f-K*U_Boundary;
    
    f_tilde(BoundaryNodes) = [];
    
    %Reduce system
    K_tilde = K;
    K_tilde(:,BoundaryNodes) = [];
    K_tilde(BoundaryNodes,:) = [];

    

    U_tilde = K_tilde\f_tilde;
    U=U_tilde;
    
    %Reinsert the conditional displacements into displacement vector.
    BoundaryNodes=sort(BoundaryNodes);  %The vector containing the boundary nodes has to be sorted for the following method

    for i = 1:length(BoundaryNodes)
       U = [U(1:(BoundaryNodes(i)-1)); U_Boundary(BoundaryNodes(i)); U((BoundaryNodes(i)):length(U))];
    end

end
