function [t, U_dyn] = DynamicCMSFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_dot_0,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);
if ~isempty(AdditionalModes)
    AdditionalModes(boundaryNodes,:)=[];
end
[M_tilde_ii,D_tilde_ii,K_tilde_ii,M_tilde_ie,D_tilde_ie,K_tilde_ie,V_cms] = CMS(K,M,D,NumberOfModes, AdditionalModes, NodeGrid);
invM_tilde_ii = inv(M_tilde_ii);
transV = V_cms';
q_0=V_cms\U_0;
q_0(1:length(boundaryNodes))=[];
q_dot_0=V_cms\U_dot_0;
q_dot_0(1:length(boundaryNodes))=[];


X_0 = [q_0; q_dot_0];

A = [zeros(size(M_tilde_ii\K_tilde_ii)) eye(size(M_tilde_ii\K_tilde_ii));
     -M_tilde_ii\K_tilde_ii -M_tilde_ii\D_tilde_ii];

t_span = [0, 10];
opt = odeset('MaxStep',1e-1);
[t, X] = ode23s(@(t,X) CMSTimeStepIntegration(t,A,X,invM_tilde_ii,M_tilde_ie,D_tilde_ie,K_tilde_ie,transV,NodeGrid),t_span,X_0,opt);


U_e = zeros(NodeGrid(end,end),length(t));
for i = 1:length(t)
    U_e(:,i) = PositionBoundaryCondition(NodeGrid,t(i));

end
U_e = U_e(boundaryNodes,:);

q_dyn = X(:,1:end/2).';

q_dyn = [U_e; q_dyn]; 
U_dyn = V_cms*q_dyn;
U_dyn(1:length(boundaryNodes),:)=[];
U_dyn = BoundaryReinsertion(NodeGrid,t,boundaryNodes,U_dyn);
end