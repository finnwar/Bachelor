function dXdt = timeStepIntegrationModal(t,A,B,C,D,E,X,NodeGrid)
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Extract only external dofs
    U_e = U_e(BoundaryNodes);
    U_e_dot = U_e_dot(BoundaryNodes);
    U_e_ddot = U_e_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_hat = B*f-C*U_e-D*U_e_dot-E*U_e_ddot;

    dXdt = A*X + [zeros(size(f_hat)); f_hat];
end