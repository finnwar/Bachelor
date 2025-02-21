function dXdt = CMSTimeStepIntegration(t,A,X,invM_tilde_ii,M_tilde_ie,D_tilde_ie,K_tilde_ie,transV,NodeGrid)
    
    % Get Boundary conditions
    [U_e, BoundaryNodes, U_e_dot, U_e_ddot] = PositionBoundaryCondition(NodeGrid,t);
    [f,~] = ForceBoundaryCondition(NodeGrid,t);
    f_e = f(BoundaryNodes);
    f(BoundaryNodes)=[];
    f_i = f;
    f=[f_e;f_i];
    f_tilde = transV*f;
    f_tilde(1:length(BoundaryNodes))=[];
    
    U_e=U_e(BoundaryNodes);
    U_e_dot=U_e_dot(BoundaryNodes);
    U_e_ddot=U_e_ddot(BoundaryNodes);
    
    f_hat = invM_tilde_ii*(f_tilde-M_tilde_ie*U_e_ddot-D_tilde_ie*U_e_dot-K_tilde_ie*U_e);
    dXdt = (A*X + [zeros(size(f_hat));f_hat]);
end