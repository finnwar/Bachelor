function dXdt = timeStepIntegration(t,A,X,M_ii_inv,M_ii_invK_ie,M_ii_invD_ie,M_ii_invM_ie,NodeGrid)
    
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_e = U_e(BoundaryNodes);
    U_e_dot = U_e_dot(BoundaryNodes);
    U_e_ddot = U_e_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_bar = M_ii_inv*f-M_ii_invK_ie*U_e-M_ii_invD_ie*U_e_dot-M_ii_invM_ie*U_e_ddot;

    dXdt = A*X + [zeros(size(f_bar)); f_bar];
end