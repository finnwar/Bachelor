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

NumberOfElementsX = 80;
NumberOfElementsY = 10;

%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 0;
beta = 0;

D = alpha*K+beta*M;

K = sparse(K);
M = sparse(M);
D = sparse(D);

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
[t_cms40, U_dyn_cms40] = DynamicCMSFEM(K,M,D,NodeGrid,40,[]);
toc;
%%
tic;
[t_cms20, U_dyn_cms20] = DynamicCMSFEM(K,M,D,NodeGrid,20,[]);
toc;

%%
tic;
[t_cms10, U_dyn_cms10] = DynamicCMSFEM(K,M,D,NodeGrid,10,[]);
toc;
%%
tic;
[t_cms4, U_dyn_cms4] = DynamicCMSFEM(K,M,D,NodeGrid,4,[]);
toc;
%%
tic;
[t_cms1, U_dyn_cms1] = DynamicCMSFEM(K,M,D,NodeGrid,1,[]);
toc;
%% Stress Calculation
tic;
[nodeStress_cms1, abscissaStress_cms1] = StressCalculation(U_dyn_cms1,t_cms1,nu,E,NodePositionTable,NodeTable);
toc;
%%
tic;
[nodeStress_cms10, abscissaStress_cms10] = StressCalculation(U_dyn_cms10,t_cms10,nu,E,NodePositionTable,NodeTable);
toc;
%%
tic;
[nodeStress_dir, abscissaStress_dir] = StressCalculation(U_dyn_dir,t_dir,nu,E,NodePositionTable,NodeTable);
toc;
%% Visualisation

PatchPlot('Direct Time Integration',U_dyn_dir,t_dir,nodeStress_dir,NodePosition,NodeTable,NumberOfElementsX,NumberOfElementsY, ...
                                            length_end, length_middle, thickness_end, thickness_middle)
%%
PatchPlot('1 Eigenmode',U_dyn_cms1,t_cms1,nodeStress_cms1,NodePosition,NodeTable,NumberOfElementsX,NumberOfElementsY, ...
                                            length_end, length_middle, thickness_end, thickness_middle)
%%
PatchPlot('10 Eigenmode',U_dyn_cms10,t_cms10,nodeStress_cms10,NodePosition,NodeTable,NumberOfElementsX,NumberOfElementsY, ...
                                            % length_end, length_middle, thickness_end, thickness_middle)
%%

ErrorPlot(U_dyn_dir,t_dir,U_dyn_cms1,t_cms1,'CMS 1')
ErrorPlot(U_dyn_dir,t_dir,U_dyn_cms4,t_cms4,'CMS 4')
ErrorPlot(U_dyn_dir,t_dir,U_dyn_cms10,t_cms10,'CMS 10')
ErrorPlot(U_dyn_dir,t_dir,U_dyn_cms20,t_cms20,'CMS 20')
ErrorPlot(U_dyn_dir,t_dir,U_dyn_cms40,t_cms40,'CMS 40')

%%
[~, MASE(1), ~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms1,U_dyn_cms1);
[~, MASE(2), ~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms4,U_dyn_cms4);
[~, MASE(3), ~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms10,U_dyn_cms10);
[~, MASE(4), ~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms20,U_dyn_cms20);
[~, MASE(5), ~,~] = ErrorCalculation(t_dir,U_dyn_dir,t_cms40,U_dyn_cms40);

figure
plot([1 4 10 20 40],MASE)
