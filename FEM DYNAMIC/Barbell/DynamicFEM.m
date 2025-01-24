function [t, U_dyn] = DynamicFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_boundary, boundaryNodes,~,~] = PositionBoundaryCondition(NodeGrid,0);
boundaryNodes = sort(boundaryNodes);

% Reconfiguration of matrices
% [ A_bb A_bf;
%   A_fb A_ff]
K_tilde = MatrixReconfiguration(K, boundaryNodes);

M_tilde = MatrixReconfiguration(M, boundaryNodes);

D_tilde = MatrixReconfiguration(D, boundaryNodes);



U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end