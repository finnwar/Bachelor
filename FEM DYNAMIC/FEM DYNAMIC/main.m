


%Simulation eines einfachen Stabes mithilfe einer 2D-FEM

[K, M, f, U] = Static_FEM(100, 5);
%%
%Eigenmdoes
[eigenVectors,eigenFrequencies] = eig(K,M,'vector');

[eigenFrequencies, index] = sort(eigenFrequencies);
eigenVectors = eigenVectors(:,index);
eigenFrequencies = sqrt(eigenFrequencies)/(2*pi);

%Damping Matrix via Rayleigh-Coefficients. PROVISIONAL STAND IN
alpha = 0.000005;
beta = 0.000025;

D = alpha * M + beta * K;

% Model Reduction
Phi = [eigenVectors(:,1:100)];

M_tilde = Phi'*M*Phi;
K_tilde = Phi'*M*Phi;
D_tilde = Phi'*D*Phi;
f_tilde = Phi'*f;

% The modal reduction produces under rank matrices when unfit modes are
% chosen. Maybe not lin. ind. from eigenmodes?
% Transient Analysis

%%
M_tilde_inv = inv(M_tilde);

U_0 = U;
q_0 = Phi\U_0;

qdot_0 = zeros(length(q_0),1);
qddot_0 = zeros(length(q_0),1);



f_tilde_time = zeros(size(f_tilde));
Y0 = [q_0;qdot_0];
A = [zeros(size(K_tilde)) eye(size(K_tilde));
     -M_tilde\K_tilde zeros(size(K_tilde))];

% Implementation of Damping Term leads to problems

% A = [zeros(size(K_tilde)) eye(size(K_tilde));
%      -M_tilde\K_tilde -M_tilde\D_tilde];
b = [zeros(size(f_tilde_time));M_tilde\f_tilde_time];
%% 
opt = odeset(MaxStep = 1e-1);
tspan = linspace(0, 60, 6000);

[t,Ysol] = ode23s(@(t,Y) timeStepIntegration(t,Y,A,b), tspan, Y0 );


U_sol = Phi* Ysol(:,1:length(Ysol(1,:))/2)';
figure
plot(t,U_sol(end,:)')

%%
% Compute the Fourier Transform of the solution
amplitude = fft(U_sol(end,:));

% Compute the frequency axis
Fs = 1 / (t(2) - t(1));  % Sampling frequency
n = length(t);           % Number of samples
frequencyRange = (0:n-1)*(Fs/n);      % Frequency range

% Plot the magnitude of the Fourier Transform
figure;
loglog(frequencyRange, abs(amplitude));
title('Fourier Transform of ODE Solution');
xlabel('Frequency (Hz)');
ylabel('Magnitude');