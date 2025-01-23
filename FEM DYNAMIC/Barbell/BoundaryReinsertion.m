function U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn)
U_dyn = zeros(NodeGrid(end,end),length(t));
%Reinsert the conditional displacements into displacement vector.
BoundaryNodes=sort(BoundaryNodes);  %The vector containing the boundary nodes has to be sorted for the following method

for j = 1:length(t)
    U_Boundary = PositionBoundaryCondition(NodeGrid,t(j));
    temp = U_tilde_dyn(:,j);
    for i = 1:length(BoundaryNodes)
        temp = [temp(1:(BoundaryNodes(i)-1)); U_Boundary(BoundaryNodes(i)); temp(BoundaryNodes(i):length(temp))];
    end
    U_dyn(:,j)=temp;
end
end