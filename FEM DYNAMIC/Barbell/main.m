%%=========================================================================
% Main Script of Bachelor Thesis of Finn Warlimont
%%=========================================================================

% Dimensions

length_total = 0.400;     %Barbell length [m]
length_middle = 0.200;    %Middle part length [m]
length_end = (length_total-length_middle)/2;

thickness = 0.01;

thickness_end = 0.2;    %Thickness of both ends [m]
thickness_middle = 0.1;   %Thickness of thin middle part [m]

% Material Constants

rho = 7850;      %Density [kg/m^3]
E = 210e9;       %Young's Modulus [N/m^2]
nu = 0.3;        %Contraction [-]

% Gravity

g = 9.81;       %Gravitational constant [m/s^2]

% Mesh-Resolution

NumberOfElementsX = 20;
NumberOfElementsY = 10;

%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 0;
beta = 1e-3;

D = alpha*K+beta*M;


[Phi_vM,PhiX,PhiY,PhiXY,PhiXnode,PhiYnode,PhiXYnode,PhiXcenter,PhiYcenter,PhiXYcenter] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable,nu,E);

%% Solve static FEM with boundary conditions
f = zeros(NodeGrid(end,end),1);
f(end) = 1000;
U_static = StaticFEM(K,f,NodeGrid);
U_mass = StaticFEM(K,g*M*ones(size(U_static)),NodeGrid);
% PlotDisplacement(U_static,NodeGrid,NodePosition)

%% Transient response
tic;
[t_dir, U_dyn_dir] = DynamicFEM(K,M,D,NodeGrid);
toc;
%%
Kernel = null(K);
%%
K = sparse(K);
M = sparse(M);
D = sparse(D);



    sigXref = PhiX*U_dyn_dir;
    sigYref = PhiY*U_dyn_dir;
    tauXYref = PhiXY*U_dyn_dir;

    sigXnoderef = PhiXnode*U_dyn_dir;
    sigYnoderef = PhiYnode*U_dyn_dir;
    tauXYnoderef = PhiXYnode*U_dyn_dir;
    
    sigXcenterref = PhiXcenter*U_dyn_dir;
    sigYcenterref = PhiYcenter*U_dyn_dir;
    tauXYcenterref = PhiXYcenter*U_dyn_dir;
%% Convergence Study of Stresses
numberOfModes = zeros(4,1);
total_error=numberOfModes;
mean_rel_error=numberOfModes;


maximum_modes=441;
interval=20;
for i = 1:interval:maximum_modes
    numberOfModes((i-1)/interval+1)=i;
    j=(i-1)/interval+1;
    [t,U_dyn_cms]=DynamicCMSFEM(K,M,D,NodeGrid,i,[]);
    [~,~,total_error((i-1)/interval+1),mean_rel_error((i-1)/interval+1)]=ErrorCalculation(t_dir,U_dyn_dir,t,U_dyn_cms);


    sigX = PhiX*U_dyn_cms;
    sigY = PhiY*U_dyn_cms;
    tauXY = PhiXY*U_dyn_cms;

    sigXnode = PhiXnode*U_dyn_cms;
    sigYnode = PhiYnode*U_dyn_cms;
    tauXYnode = PhiXYnode*U_dyn_cms;

    sigXcenter = PhiXcenter*U_dyn_cms;
    sigYcenter = PhiYcenter*U_dyn_cms;
    tauXYcenter = PhiXYcenter*U_dyn_cms;

    [~, ~, absStressErrorX(j), relStressErrorX(j),~]=StressErrorCalculation(t_dir,sigXref,t,sigX);
    [~, ~, absStressErrorY(j), relStressErrorY(j),~]=StressErrorCalculation(t_dir,sigYref,t,sigY);
    [~, ~, absStressErrorXY(j), relStressErrorXY(j),~]=StressErrorCalculation(t_dir,tauXYref,t,tauXY);

    [~, ~, absStressErrorXnode(j), relStressErrorXnode(j),~]=StressErrorCalculation(t_dir,sigXnoderef,t,sigXnode);
    [~, ~, absStressErrorYnode(j), relStressErrorYnode(j),~]=StressErrorCalculation(t_dir,sigYnoderef,t,sigYnode);
    [~, ~, absStressErrorXYnode(j), relStressErrorXYnode(j),~]=StressErrorCalculation(t_dir,tauXYnoderef,t,tauXYnode);

    [~, ~, absStressErrorXcenter(j), relStressErrorXcenter(j),~]=StressErrorCalculation(t_dir,sigXcenterref,t,sigXcenter);
    [~, ~, absStressErrorYcenter(j), relStressErrorYcenter(j),~]=StressErrorCalculation(t_dir,sigYcenterref,t,sigYcenter);
    [~, ~, absStressErrorXYcenter(j), relStressErrorXYcenter(j),~]=StressErrorCalculation(t_dir,tauXYcenterref,t,tauXYcenter);


end


%% Plotting
figure('Name','Convergence of Relative Error of Stress at Element Center')

plot(numberOfModes,relStressErrorXcenter,'Color','blue')
hold on;
plot(numberOfModes,relStressErrorYcenter,'Color','red')
plot(numberOfModes,relStressErrorXYcenter,'Color','black')
xlabel('Number of Modes')
ylabel('Relative Error')
title('Convergence of Stress Evaluated at Element Center')
legend('Normal Stress X','Normal Stress Y','Shear Stress')


figure('Name','Convergence of Relative Error of Stress at Gauss Points')

plot(numberOfModes,relStressErrorX,'Color','blue')
hold on;
plot(numberOfModes,relStressErrorY,'Color','red')
plot(numberOfModes,relStressErrorXY,'Color','black')
xlabel('Number of Modes')
ylabel('Relative Error')
title('Convergence of Stress Evaluated at Gauss Points')
legend('Normal Stress X','Normal Stress Y','Shear Stress')


figure('Name','Convergence of Relative Error of Stress at Nodes')

plot(numberOfModes,relStressErrorXnode,'Color','blue')
hold on;
plot(numberOfModes,relStressErrorYnode,'Color','red')
plot(numberOfModes,relStressErrorXYnode,'Color','black')
xlabel('Number of Modes')
ylabel('Relative Error')
title('Convergence of Stress Evaluated at Nodes')
legend('Normal Stress X','Normal Stress Y','Shear Stress')
