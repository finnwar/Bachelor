


%Simulation eines einfachen Stabes mithilfe einer 2D-FEM

[K, M, f, U] = Static_FEM;

%Eigenmdoes
[eigenVectors,eigenFrequencies] = eig(K,M,'vector');

[eigenFrequencies, index] = sort(eigenFrequencies);
eigenVectors = eigenVectors(:,index);


%Damping Matrix via Rayleigh-Coefficients
alpha = 0.005;
beta = 5;

D = alpha * M + beta * K;

% Model Reduction
Phi = [eigenVectors(:,1:10)];

M_tilde = Phi'*M*Phi;
K_tilde = Phi'*M*Phi;
D_tilde = Phi'*D*Phi;
f_tilde = Phi'*f;

% The modal redcution produces under rank matrices when unfit modes are
% chosen. Maybe not lin. ind. from eigenmodes?
% Transient Analysis

%%
M_tilde_inv = inv(M_tilde);
%opt = odeset('maxstep',1e-9);
tspan = [0 60];

U_0 = U;
q_0 = Phi\U_0;

qdot_0 = zeros(length(q_0),1);
qddot_0 = zeros(length(q_0),1);



f_tilde_time = zeros(size(f_tilde));
Y0 = [q_0;qdot_0];
A = [zeros(size(K_tilde)) eye(size(K_tilde));
     -M_tilde\K_tilde zeros(size(K_tilde))];            % Implementation of Damping Term leads to problems
b = [zeros(size(f_tilde_time));M_tilde\f_tilde_time];
%% 

[t,Ysol] = ode23t(@(t,Y) timeStepIntegration(t,Y,A,b), tspan, Y0);


U_sol = Phi* Ysol(:,1:length(Ysol(1,:))/2)';
figure
plot(t,U_sol(end/4,:)')