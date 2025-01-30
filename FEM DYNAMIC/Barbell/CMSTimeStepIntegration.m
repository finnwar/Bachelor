function dXdt = CMSTimeStepIntegration(t,A,X,transposedPhi,invM_ii,B_M,B_D,B_K,NodeGrid)
    
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_e = U_e(BoundaryNodes);
    U_e_dot = U_e_dot(BoundaryNodes);
    U_e_ddot = U_e_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_tilde = transposedPhi*(f-B_M*U_e_ddot-B_D*U_e_dot-B_K*U_e);
    
    dXdt = invM_ii*(A*X + [zeros(size(f_tilde));f_tilde]);
end