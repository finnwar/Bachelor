function [t, U_dyn] = DynamicModalFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_0_dot,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);

[Omega, Phi, D_tilde, D_bar, M_bar, invK_FFK_FB] = ModalReduction(K,M,D,NumberOfModes,AdditionalModes,NodeGrid);

U_F0 = U_0;
U_F0(boundaryNodes) = [];
U_B0 = U_0(boundaryNodes);

U_F0_dot = U_0_dot;
U_F0_dot(boundaryNodes) = [];
U_B0_dot = U_0_dot(boundaryNodes);



q_0 = Phi\(U_F0+invK_FFK_FB*U_B0);
q_0_dot = Phi\(U_F0_dot+invK_FFK_FB*U_B0_dot);

A = [zeros(size(Omega)) eye(size(Omega));
     -Omega -D_tilde]; 
transposedPhi = Phi.';


t_span = [0 10];
[t, X] = ode23(@(t,X) ModalTimeStepIntegration(t,A,X,transposedPhi,M_bar,D_bar,NodeGrid),t_span,[q_0;q_0_dot]);


U_dyn = BoundaryReinsertion(NodeGrid, t, boundaryNodes, Phi*X(:,1:(end/2)).');



end