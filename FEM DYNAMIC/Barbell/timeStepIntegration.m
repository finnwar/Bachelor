function dXdt = timeStepIntegration(t,A,X,M_tilde_inv,NodeGrid, K,M,D,transposedPhi)
    
    % Get Boundary conditions
    [U_b, BoundaryNodes, U_b_dot, U_b_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_b(BoundaryNodes) = [];
    U_b_dot(BoundaryNodes) = [];
    U_b_ddot(BoundaryNodes) = [];
    f(BoundaryNodes) = [];

    f_tilde = f-K*U_b-D*U_b_dot-M*U_b_ddot;
    
    dXdt = A*X + [zeros(size(f_tilde)); M_tilde_inv*f_tilde];
end