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

%% Solving with Modal Reduction
NumberOfModes = NodeGrid(end,end)/10;


[K_hat, M_hat, D_hat, ~,Phi, ~] = ModalReduction(K_tilde, M_tilde, D_tilde, f_tilde, NumberOfModes);
M_hat_inv = inv(M_hat);
A = [zeros(size(K_hat)) eye(size(K_hat));
     -M_hat\K_hat -M_hat\D_hat];
%Initial Conditions

[u0,~,u_dot0] = PositionBoundaryCondition(NodeGrid,0);
u0(BoundaryNodes) = [];
u_dot0(BoundaryNodes) = [];

q0 = Phi\u0;
q_dot0 = Phi\u_dot0;

X0 = [q0;q_dot0];

tspan = [0 1];
opt = odeset('Maxstep',1e-3);
[t,Xsol] = ode15s(@(t,X) timeStepIntegration(t,A,X,M_hat_inv,NodeGrid,K_tilde,M_tilde,D_tilde,Phi.'), tspan, X0,opt);




% %% Solving without Modal Reduction
% A = [zeros(size(K_tilde)) eye(size(K_tilde));
%      -M_tilde\K_tilde -M_tilde\D_tilde];            
% 
% [u0,~,u_dot0] = PositionBoundaryCondition(NodeGrid,0);
% 
% u0(BoundaryNodes) = [];
% u_dot0(BoundaryNodes) = [];
% 
% %Initial Conditions
% X0 = [u0;u_dot0];
% 
% A = [zeros(size(K_tilde)) eye(size(K_tilde)); -M_tilde\K_tilde -M_tilde\D_tilde];
% 
% %% ODE Solver
% tspan = [0 1];
% [t,Xsol] = ode15s(@(t,X) timeStepIntegration(t,A,X,K_tilde,M_tilde,invM_tilde,D_tilde,NodeGrid), tspan, X0);

%% Post-Processing

%No modal reduction:
% U_tilde_dyn = Xsol(:,1:(end/2)).';
%With modal reduction:
U_tilde_dyn = Phi*Xsol(:,1:(end/2)).';



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