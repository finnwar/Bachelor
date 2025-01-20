function dXdt = timeStepIntegration(t,A,X,M_tilde_inv,NodeGrid, K,M,D,transposedPhi)
    
    % Get Boundary conditions
    [U, BoundaryNodes, U_dot, U_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    
    % Reduce to eliminate bound degrees of freedom
    U(BoundaryNodes) = [];
    U_dot(BoundaryNodes) = [];
    U_ddot(BoundaryNodes) = [];
    f(BoundaryNodes) = [];

    f_tilde = f-K*U-D*U_dot-M*U_ddot;
    if nargin == 9
        f_hat = transposedPhi*f_tilde;
    else
        f_hat = f_tilde;
    end

    
    dXdt = A*X + [zeros(size(f_hat)); M_tilde_inv*f_hat];
end