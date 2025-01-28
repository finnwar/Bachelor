function dXdt = ModalTimeStepIntegration(t,A,X,transposedPhi,invM_tilde_ii,M_tilde_ie,D_tilde_ie,K_tilde_ie,NodeGrid)
    
    % Get Boundary conditions
    [U_b, BoundaryNodes, U_b_dot, U_b_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_b = U_b(BoundaryNodes);
    U_b_dot = U_b_dot(BoundaryNodes);
    U_b_ddot = U_b_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_tilde = transposedPhi*f - M_tilde_ie*U_b_ddot - D_tilde_ie*U_b_dot - K_tilde_ie*U_b;
    
    dXdt = A*X + [zeros(size(f_tilde));invM_tilde_ii*f_tilde];
end