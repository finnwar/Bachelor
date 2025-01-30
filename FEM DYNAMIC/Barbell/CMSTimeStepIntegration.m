function dXdt = CMSTimeStepIntegration(t,A,X,invM_tildetransposedPhi,invM_tildeM_tilde_ie,invM_tildeD_tilde_ie, invM_tildeK_tilde_ie,NodeGrid)
    
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_e = U_e(BoundaryNodes);
    U_e_dot = U_e_dot(BoundaryNodes);
    U_e_ddot = U_e_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_tilde = invM_tildetransposedPhi*f - invM_tildeM_tilde_ie*U_e_ddot - invM_tildeD_tilde_ie*U_e_dot - invM_tildeK_tilde_ie*U_e;
    
    dXdt = A*X + [zeros(size(f_tilde));f_tilde];
end