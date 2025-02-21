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

NumberOfElementsX = 30;
NumberOfElementsY = 10;

%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 1e-10;
beta = 1e-4;

D = alpha*K+beta*M;


[Phi_vM,PhiX,PhiY,PhiXY] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable,nu,E);
%% Solve static FEM with boundary conditions
f = zeros(NodeGrid(end,end),1);
f(end) = -1000;
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
[t_cms100, U_dyn_cms100] = DynamicCMSFEM(K,M,D,NodeGrid,100,[]);
toc;

%%
tic;
[t_cms270, U_dyn_cms270] = DynamicCMSFEM(K,M,D,NodeGrid,270,[]);
toc;
%% Visualisation 
%% Nodal Approach
PatchPlot('Nodal Approach',U_dyn_dir,t_dir,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('CMS 5',U_dyn_cms5,t_cms5,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 
%% 

PatchPlot('CMS 20',U_dyn_cms20,t_cms20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('Modal 20',U_dyn_mod20,t_mod20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

%% 

PatchPlot('Modal 100',U_dyn_mod100,t_mod100,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('CMS 100',U_dyn_cms100,t_cms100,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('Modal 270',U_dyn_mod270,t_mod270,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('CMS 270',U_dyn_cms270,t_cms270,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% 

PatchPlot('Static',U_static,[1],Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)
%% Error Calculation

% Deformation Error

[~,~,uError(1),~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms5,U_dyn_cms5);
[~,~,uError(2),~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms20,U_dyn_cms20);
[~,~,uError(3),~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms100,U_dyn_cms100);
[~,~,uError(4),~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms270,U_dyn_cms270);

% %% 
% [~,~,sError(1),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_cms5,U_dyn_cms5);
% [~,~,sError(2),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_mod5,U_dyn_mod5);
% [~,~,sError(3),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_cms20,U_dyn_cms20);
% [~,~,sError(4),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_mod20,U_dyn_mod20);
% [~,~,sError(5),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_cms100,U_dyn_cms100);
% [~,~,sError(6),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_mod100,U_dyn_mod100);
% [~,~,sError(7),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_cms270,U_dyn_cms270);
% [~,~,sError(8),~] = StressErrorCalculation(t_dir,U_dyn_dir,t_mod270,U_dyn_mod270);








