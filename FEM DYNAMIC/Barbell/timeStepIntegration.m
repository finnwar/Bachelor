function dXdt = timeStepIntegration(t,A,X,K,M,D,NodeGrid)
    
    % Get Boundary conditions
    [U, BoundaryNodes, U_dot, U_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U(BoundaryNodes) = [];
    U_dot(BoundaryNodes) = [];
    U_ddot(BoundaryNodes) = [];
    f(BoundaryNodes) = [];

    f_tilde = f-K*U-D*U_dot-M*U_ddot;
    
    
    

    dXdt = A*X + [zeros(size(f_tilde)); f_tilde];
end