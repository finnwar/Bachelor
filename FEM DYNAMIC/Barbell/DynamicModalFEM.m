function [t, U_dyn] = DynamicModalFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_0_dot,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);

[V_cms,Phi, M_tilde,D_tilde,K_tilde] = ModalReduction(K,M,D,NumberOfModes,AdditionalModes,NodeGrid);


[M_tilde, M_tilde_ee, M_tilde_ei, M_tilde_ie, M_tilde_ii]= MatrixReconfiguration(M, boundaryNodes);
[D_tilde, D_tilde_ee, D_tilde_ei, D_tilde_ie, D_tilde_ii]= MatrixReconfiguration(M, boundaryNodes);
[K_tilde, K_tilde_ee, K_tilde_ei, K_tilde_ie, K_tilde_ii]= MatrixReconfiguration(M, boundaryNodes);

A = [zeros(size(K_tilde_ii)) eye(size(K_tilde_ii));
     -M_tilde_ii\K_tilde_ii -M_tilde_ii\D_tilde_ii];

q_0 = V_cms\U_0;
q_0(boundaryNodes) = [];

q_0_dot = V_cms\U_0_dot;
q_0_dot(boundaryNodes) = [];


t_span = [0 10];
[t, X] = ode23(@(t,X) ModalTimeStepIntegration(t,A,X,Phi.',M_tilde_ie,D_tilde_ie,K_tilde_ie,NodeGrid),t_span,[q_0;q_0_dot]);


U_dyn = BoundaryReinsertion(NodeGrid, t, boundaryNodes, Phi*X(:,1:(end/2)).');



end