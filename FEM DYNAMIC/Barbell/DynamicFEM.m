function [t, U_dyn] = DynamicFEM(K,M,D,NodeGrid,NumberOfModes,AdditionalModes)

[U_Boundary, BoundaryNodes] = PositionBoundaryCondition(NodeGrid,0);
[f,~] = ForceBoundaryCondition(NodeGrid,0);

f_tilde = f-K*U_Boundary;

f_tilde(BoundaryNodes) = [];

% Eliminate constraints
K_tilde = K;
K_tilde(:,BoundaryNodes) = [];
K_tilde(BoundaryNodes,:) = [];

D_tilde = D;
D_tilde(:,BoundaryNodes) = [];
D_tilde(BoundaryNodes,:) = [];

M_tilde = M;
M_tilde(:,BoundaryNodes) = [];
M_tilde(BoundaryNodes,:) = [];

M_tilde_inv = inv(M_tilde);

%% Solving with Modal Reduction


if nargin > 4
    if nargin == 6
        AdditionalModes(BoundaryNodes,:) = [];
        [K_hat, M_hat, D_hat, ~,Phi, ~] = ModalReduction(K_tilde, M_tilde, D_tilde, f_tilde, NumberOfModes, AdditionalModes);
    elseif nargin == 5
        [K_hat, M_hat, D_hat, ~,Phi, ~] = ModalReduction(K_tilde, M_tilde, D_tilde, f_tilde, NumberOfModes);
    end
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


    %With modal reduction:
    U_tilde_dyn = Phi*Xsol(:,1:(end/2)).';


else
    %% Solving without Modal Reduction

    [u0,~,u_dot0] = PositionBoundaryCondition(NodeGrid,0);

    u0(BoundaryNodes) = [];
    u_dot0(BoundaryNodes) = [];

    %Initial Conditions
    X0 = [u0;u_dot0];

    A = [zeros(size(K_tilde)) eye(size(K_tilde));
        -M_tilde\K_tilde -M_tilde\D_tilde];

    %% ODE Solver
    tspan = [0 10];
    [t,Xsol] = ode15s(@(t,X) timeStepIntegration(t,A,X,M_tilde_inv,NodeGrid,K_tilde,M_tilde,D_tilde), tspan, X0);

    % Post-Processing

    % No modal reduction:
    U_tilde_dyn = Xsol(:,1:(end/2)).';

end

U_dyn = BoundaryReinsertion(NodeGrid, t, BoundaryNodes, U_tilde_dyn);
% U_dyn_dot = BoundaryReinsertion(NodeGrid,t,BoundaryNodes,Xsol(:,(end/2+1):end),)


end