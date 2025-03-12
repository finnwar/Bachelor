function [t, U_dyn] = DynamicFEM(K,M,D,NodeGrid)

[U_Boundary, BoundaryNodes] = PositionBoundaryCondition(NodeGrid,0);
[f,~] = ForceBoundaryCondition(NodeGrid,0);

[M,M_ee, M_ei,M_ie, M_ii] = MatrixReconfiguration(M,BoundaryNodes);

[K,K_ee, K_ei,K_ie, K_ii] = MatrixReconfiguration(K,BoundaryNodes);

[D,D_ee, D_ei,D_ie, D_ii] = MatrixReconfiguration(D,BoundaryNodes);

M_ii_inv = inv(M_ii);
M_ii_invK_ie = M_ii_inv*K_ie;
M_ii_invD_ie = M_ii_inv*D_ie;
M_ii_invM_ie = M_ii_inv*M_ie;
    %% Solving without Modal Reduction

    [u_i0,~,u_i_dot0] = PositionBoundaryCondition(NodeGrid,0);

    u_i0(BoundaryNodes) = [];
    u_i_dot0(BoundaryNodes) = [];

    %Initial Conditions
    X0 = [u_i0;u_i_dot0];

    A = [zeros(size(K_ii)) eye(size(K_ii));
        -M_ii\K_ii -M_ii\D_ii];

    %% ODE Solver
    opt = odeset('MaxStep',1e-1);
    tspan = linspace(0,10,1000);
    [t,Xsol] = ode15s(@(t,X) timeStepIntegration(t,A,X,M_ii_inv,M_ii_invK_ie,M_ii_invD_ie,M_ii_invM_ie,NodeGrid), tspan, X0,opt);

    % Post-Processing

    % No modal reduction:
    U_tilde_dyn = Xsol(:,1:(end/2)).';

U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end