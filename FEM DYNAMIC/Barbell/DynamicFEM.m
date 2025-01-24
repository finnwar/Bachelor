function [t, U_dyn] = DynamicModalFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_boundary, boundaryNodes,~,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);



U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end