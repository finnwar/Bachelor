%%=========================================================================
% Main Script of Bachelor Thesis of Finn Warlimont
%%=========================================================================

% Dimensions

length_total = 0.400;     %Barbell length [m]
length_middle = 0.200;    %Middle part length [m]
length_end = (length_total-length_middle)/2;

thickness = 0.1;

thickness_end = 0.2;    %Thickness of both ends [m]
thickness_middle = 0.1;   %Thickness of thin middle part [m]

% Material Constants

rho = 7850;      %Density [kg/m^3]
E = 210e9;       %Young's Modulus [N/m^2]
nu = 0.3;        %Contraction [-]

% Gravity

g = 9.81;       %Gravitational constant [m/s^2]

% Mesh-Resolution

NumberOfElementsX = 1;
NumberOfElementsY = 1;



%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 0;
beta = 0;


D = alpha*K+beta*M;
%% Solve static FEM with boundary conditions

U_static = StaticFEM(K,ForceBoundaryCondition(NodeGrid,0),NodeGrid);
U_mass = StaticFEM(K,g*M*ones(size(U_static)),NodeGrid);
% PlotDisplacement(U_static,NodeGrid,NodePosition)

%% Transient response
[t_modal10, U_dyn_modal10] = DynamicFEM(K,M,D,NodeGrid,10);
disp('Modal10')
% [t_modal40, U_dyn_modal40] = DynamicFEM(K,M,D,NodeGrid,10);
% disp('Modal40')
% [t_modal40add, U_dyn_modal40add] = DynamicFEM(K,M,D,NodeGrid,40, [U_static/norm(U_static), U_mass/norm(U_mass)]);
% disp('Modal40add')
% [t_modal10add, U_dyn_modal10add] = DynamicFEM(K,M,D,NodeGrid,10, [U_static/norm(U_static), U_mass/norm(U_mass)]);
% disp('Modal10add')
% [t_modal10Ker, U_dyn_modal10Ker] = DynamicFEM(K,M,D,NodeGrid,10, null(K));
% disp('Modal10Ker')
% [t_modal40Ker, U_dyn_modal40Ker] = DynamicFEM(K,M,D,NodeGrid,40, null(K));
% disp('Modal40Ker')
% 
[~,posBoundNodes,~] = PositionBoundaryCondition(NodeGrid,0);

[t_modalFull, U_dyn_modalFull] = DynamicFEM(K,M,D,NodeGrid,NodeGrid(end,end)-length(posBoundNodes));
% disp('ModalFull')
[t_direct, U_dyn_direct] = DynamicFEM(K,M,D,NodeGrid);
disp('Direct')





