function [t, U_dyn] = DynamicCMSFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_0_dot,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);
AdditionalModes(boundaryNodes,:) = [];
[V_cms,Phi, M_tilde,D_tilde,K_tilde,K_ii_invK_ie] = CMS(K,M,D,NumberOfModes,AdditionalModes,NodeGrid);


[M_tilde, M_tilde_ee, M_tilde_ei, M_tilde_ie, M_tilde_ii]= MatrixReconfiguration(M_tilde, 1:length(boundaryNodes));
[D_tilde, D_tilde_ee, D_tilde_ei, D_tilde_ie, D_tilde_ii]= MatrixReconfiguration(D_tilde, 1:length(boundaryNodes));
[K_tilde, K_tilde_ee, K_tilde_ei, K_tilde_ie, K_tilde_ii]= MatrixReconfiguration(K_tilde, 1:length(boundaryNodes));

invM_tilde_ii = inv(M_tilde_ii);
invM_tildetransposedPhi = invM_tilde_ii*Phi.';
invM_tildeM_tilde_ie = invM_tilde_ii*M_tilde_ie;
invM_tildeD_tilde_ie = invM_tilde_ii*D_tilde_ie;
invM_tildeK_tilde_ie = invM_tilde_ii*K_tilde_ie;

A = [zeros(size(K_tilde_ii)) eye(size(K_tilde_ii));
     -M_tilde_ii\K_tilde_ii -M_tilde_ii\D_tilde_ii];



U_0 = [U_0(boundaryNodes); U_0(setdiff(1:length(U_0),boundaryNodes))];
U_0_dot = [U_0_dot(boundaryNodes); U_0_dot(setdiff(1:length(U_0),boundaryNodes))];


q_0 = V_cms\U_0;
q_0(1:length(boundaryNodes)) = [];


q_0_dot = V_cms\U_0_dot;
q_0_dot(1:length(boundaryNodes)) = [];


t_span = [0 10];
opt = odeset('MaxStep',1e-1);
[t, X] = ode15s(@(t,X) CMSTimeStepIntegration(t,A,X,invM_tildetransposedPhi,invM_tildeM_tilde_ie,invM_tildeD_tilde_ie,invM_tildeK_tilde_ie,NodeGrid),t_span,[q_0;q_0_dot],opt);



for i = 1:length(t)
    U_e(:,i) = PositionBoundaryCondition(NodeGrid,t(i));

end
    U_e(setdiff(1:length(U_e(:,1)),boundaryNodes),:)=[];

U_dyn = BoundaryReinsertion(NodeGrid,t,boundaryNodes,[-K_ii_invK_ie Phi]*[U_e;X(:,1:end/2).']);


end