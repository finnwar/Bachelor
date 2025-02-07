function [t, U_dyn] = DynamicCMSFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_0, boundaryNodes,U_dot_0,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);
if ~isempty(AdditionalModes)
    AdditionalModes(boundaryNodes,:)=[];
end
[invMiiPhiT,invMiiMie,invMiiDie,invMiiKie,invMiiDii,invMiiKii,V_cms,Phi] = CMS(K,M,D,NumberOfModes, AdditionalModes, NodeGrid);


q_0=V_cms\U_0;
q_0(1:length(boundaryNodes))=[];
q_dot_0=V_cms\U_dot_0;
q_dot_0(1:length(boundaryNodes))=[];


X_0 = [q_0; q_dot_0];

A = [zeros(size(invMiiKii)) eye(size(invMiiKii));
     -invMiiKii -invMiiDii];

t_span = [0 10];
opt = odeset('MaxStep',1e-1);
[t, X] = ode15s(@(t,X) CMSTimeStepIntegration(t,A,X,invMiiPhiT,invMiiMie,invMiiDie,invMiiKie,NodeGrid),t_span,X_0,opt);


U_e = zeros(NodeGrid(end,end),length(t));
for i = 1:length(t)
    U_e(:,i) = PositionBoundaryCondition(NodeGrid,t(i));

end


q_dyn = X(:,1:end/2).';

U_e = U_e(boundaryNodes,:);

U_dyn = V_cms*[U_e; q_dyn];


end