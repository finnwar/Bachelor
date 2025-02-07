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

NumberOfElementsX = 40;
NumberOfElementsY = 10;

%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 0;
beta = 0;

D = alpha*K+beta*M;

% K = sparse(K);
% M = sparse(M);
% D = sparse(D);

[Phi_vM,PhiX,PhiY,PhiXY] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable,nu,E);
%% Solve static FEM with boundary conditions

U_static = StaticFEM(K,ForceBoundaryCondition(NodeGrid,0),NodeGrid);
U_mass = StaticFEM(K,g*M*ones(size(U_static)),NodeGrid);
% PlotDisplacement(U_static,NodeGrid,NodePosition)

%% Transient response
tic;
[t_dir, U_dyn_dir] = DynamicFEM(K,M,D,NodeGrid);
toc;

%%
tic;
[t_cms20, U_dyn_cms20] = DynamicCMSFEM(K,M,D,NodeGrid,20,[]);
toc;
%%
tic;
[t_cms100, U_dyn_cms100] = DynamicCMSFEM(K,M,D,NodeGrid,100,[]);
toc;
%%
Ker = null(K);
tic;
[t_cms20K, U_dyn_cms20K] = DynamicCMSFEM(K,M,D,NodeGrid,20,Ker);
toc;
%% Visualisation 
%% Nodal Approach
ElementStressDir = StressCalculation(U_dyn_dir,t_dir,nu,E,NodePositionTable,NodeTable);
PatchPlot('Nodal Approach',U_dyn_dir,t_dir,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

%% 
ElementStressCMS20 = StressCalculation(U_dyn_cms20,t_cms20,nu,E,NodePositionTable,NodeTable);
PatchPlot('CMS 20',U_dyn_cms20,t_cms20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)