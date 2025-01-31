function [t, U_dyn] = DynamicCMSFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_0_dot,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);


[M_ii,D_ii,K_ii,K_ie,B_M,B_D,B_K,Phi,Omega] = CMS(K,M,D,NumberOfModes, AdditionalModes, NodeGrid);


U_0(boundaryNodes) = [];
U_0_dot(boundaryNodes) = [];
q_0 = Phi\U_0;
q_0_dot = Phi\U_0_dot;

X_0 = [q_0; q_0_dot];
invM_tilde_ii = inv(Phi.'*M_ii*Phi);
A = [zeros(length(Phi(1,:))) eye(length(Phi(1,:)));
     -invM_tilde_ii*Phi.'*K_ii*Phi -invM_tilde_ii*Phi.'*D_ii*Phi];

t_span = [0 10];
opt = odeset('MaxStep',1e-1);
[t, X] = ode15s(@(t,X) CMSTimeStepIntegration(t,A,X,Phi.',invM_tilde_ii,B_M,B_D,B_K,NodeGrid),t_span,X_0,opt);


U_e = zeros(NodeGrid(end,end),length(t));
for i = 1:length(t)
    U_e(:,i) = PositionBoundaryCondition(NodeGrid,t(i));

end
U_e = U_e(boundaryNodes,:);
K_ii_invK_ie = K_ii\K_ie;
U_dyn = BoundaryReinsertion(NodeGrid,t,boundaryNodes,[-K_ii_invK_ie Phi]*[U_e;X(:,1:end/2).']);


end