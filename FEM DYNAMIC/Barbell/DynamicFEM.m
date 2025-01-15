function [t, U_dyn] = DynamicFEM(K,M,D,NodeGrid)

    [U_Boundary, BoundaryNodes] = PositionBoundaryCondition(NodeGrid,0);
    [f,~] = ForceBoundaryCondition(NodeGrid,0);
    
    f_tilde = f-K*U_Boundary;
    
    f_tilde(BoundaryNodes) = [];
    
    % Eleminate constraints
    K_tilde = K;
    K_tilde(:,BoundaryNodes) = [];
    K_tilde(BoundaryNodes,:) = [];

    D_tilde = D;
    D_tilde(:,BoundaryNodes) = [];
    D_tilde(BoundaryNodes,:) = [];

    M_tilde = M;
    M_tilde(:,BoundaryNodes) = [];
    M_tilde(BoundaryNodes,:) = [];

    invM_tilde = inv(M_tilde);
tspan = [0 1];

A = [zeros(size(K_tilde)) eye(size(K_tilde));
     -M_tilde\K_tilde -M_tilde\D_tilde];            

[u0,~,u_dot0] = PositionBoundaryCondition(NodeGrid,0);

u0(BoundaryNodes) = [];
u_dot0(BoundaryNodes) = [];

X0 = [u0;u_dot0];
A = [zeros(size(K_tilde)) eye(size(K_tilde)); -M_tilde\K_tilde -M_tilde\D_tilde];
[t,Xsol] = ode15s(@(t,X) timeStepIntegration(t,A,X,K_tilde,M_tilde,invM_tilde,D_tilde,NodeGrid), tspan, X0);

U_tilde_dyn = Xsol(:,1:(end/2)).';
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