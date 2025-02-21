function dXdt = CMSTimeStepIntegration(t,A,X,invMiiPhiT,invMiiMie,invMiiDie,invMiiKie,NodeGrid)
    
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    U_e=U_e(BoundaryNodes);
    U_e_dot=U_e_dot(BoundaryNodes);
    U_e_ddot=U_e_ddot(BoundaryNodes);
    f(BoundaryNodes)=[];
    % Reduce to eliminate bound degrees of freedom
    f_hat= invMiiPhiT*f - invMiiMie * U_e_ddot - invMiiDie*U_e_dot - invMiiKie*U_e;
    dXdt = (A*X + [zeros(size(f_hat));f_hat]);
end