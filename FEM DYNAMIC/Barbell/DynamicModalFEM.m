function [t, U_dyn] = DynamicModalFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_Boundary, BoundaryNodes] = PositionBoundaryCondition(NodeGrid,0);
[f,~] = ForceBoundaryCondition(NodeGrid,0);

[M,M_ee, M_ei,M_ie, M_ii] = MatrixReconfiguration(M,BoundaryNodes);

[K,K_ee, K_ei,K_ie, K_ii] = MatrixReconfiguration(K,BoundaryNodes);

[D,D_ee, D_ei,D_ie, D_ii] = MatrixReconfiguration(D,BoundaryNodes);

    [u_i0,~,u_i_dot0] = PositionBoundaryCondition(NodeGrid,0);

    u_i0(BoundaryNodes) = [];
    u_i_dot0(BoundaryNodes) = [];

    if  ~isempty(AdditionalModes)
        AdditionalModes(BoundaryNodes,:) = [];
    end
    [Phi, Omega] = eigs(K_ii,M_ii,NumberOfModes,'smallestabs');    
    Phi = [Phi, AdditionalModes];
    K_ii_tilde = Phi.'*K_ii*Phi;
    D_ii_tilde = Phi.'*D_ii*Phi;
    M_ii_tilde = Phi.'*M_ii*Phi;
    invM_ii_tilde = inv(M_ii_tilde);

    q_i0 = Phi\u_i0;
    q_i_dot0 = Phi\u_i_dot0;

    X0 = [q_i0;q_i_dot0];
    A = [zeros(size(K_ii_tilde)) eye(size(K_ii_tilde));
         -invM_ii_tilde*K_ii_tilde -invM_ii_tilde*D_ii_tilde];
    
    B = invM_ii_tilde*Phi.';
    C = invM_ii_tilde*Phi.'*M_ie;
    D = invM_ii_tilde*Phi.'*D_ie;
    E = invM_ii_tilde*Phi.'*M_ie;
    opt = odeset('MaxStep',1e-1);
    tspan = [0, 10];
    [t,Xsol] = ode15s(@(t,X) timeStepIntegrationModal(t,A,B,C,D,E,X,NodeGrid), tspan, X0,opt);

    % Post-Processing

    % No modal reduction:
    q_dyn = Xsol(:,1:(end/2)).';
    U_tilde_dyn = Phi*q_dyn;
U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end