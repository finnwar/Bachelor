function dXdt = ModalTimeStepIntegration(t,A,X,transposedPhi,M_bar,D_bar,NodeGrid)
    
    % Get Boundary conditions
    [U_b, BoundaryNodes, U_b_dot, U_b_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U_b = U_b(BoundaryNodes);
    U_b_dot = U_b_dot(BoundaryNodes);
    U_b_ddot = U_b_ddot(BoundaryNodes);
    f(BoundaryNodes) = [];

    f_tilde = transposedPhi*(f+M_bar*U_b_ddot+D_bar*U_b_dot);
    if nargin == 9
        f_tilde = transposedPhi*f_tilde;
    end

    
    dXdt = A*X + [zeros(size(f_tilde));f_tilde];
end