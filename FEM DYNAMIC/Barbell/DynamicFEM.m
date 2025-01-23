function [t, U_dyn] = DynamicFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_boundary, boundaryNodes,~,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);

M_tilde = zeros(size(M));

M_tilde(1:length(boundaryNodes),:) = M(boundaryNodes,:);
M_tilde(:,1:length(boundaryNodes)) = M(:,boundaryNodes);
M_tilde((length(boundaryNodes)+1):end,:) = M(setdiff(1:length(M(1,:)),boundaryNodes),:);
M_tilde(:,(length(boundaryNodes)+1):end) = M(:,setdiff(1:length(M(1,:)),boundaryNodes));





U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end