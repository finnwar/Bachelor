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
NumberOfElementsY = 5;

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
[t_mod20, U_dyn_mod20] = DynamicModalFEM(K,M,D,NodeGrid,20,[]);
toc;
%%
tic;
[t_cms20, U_dyn_cms20] = DynamicCMSFEM(K,M,D,NodeGrid,20,[]);
toc;
%%
tic;
[t_cms5, U_dyn_cms5] = DynamicCMSFEM(K,M,D,NodeGrid,5,[]);
toc;
%%
tic;
[t_mod5, U_dyn_mod5] = DynamicModalFEM(K,M,D,NodeGrid,5,[]);
toc;
%%
tic;
[t_cms240, U_dyn_cms240] = DynamicCMSFEM(K,M,D,NodeGrid,240,[]);
toc;
%%
tic;
[t_mod240, U_dyn_mod240] = DynamicModalFEM(K,M,D,NodeGrid,240,[]);
toc;

%% Visualisation 
%% Nodal Approach
PatchPlot('Nodal Approach',U_dyn_dir,t_dir,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

%% 

PatchPlot('CMS 20',U_dyn_cms20,t_cms20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('Modal 20',U_dyn_mod20,t_mod20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

%% 

PatchPlot('Modal 240',U_dyn_mod240,t_mod240,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('CMS 240',U_dyn_cms240,t_cms240,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%%
[~,MASE(1),~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_mod5,U_dyn_mod5);
[~,MASE(2),~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms5,U_dyn_cms5);
[~,MASE(3),~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_mod20,U_dyn_mod20);
[~,MASE(4),~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms20,U_dyn_cms20);