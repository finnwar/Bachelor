%%=========================================================================
% Main Script of Bachelor Thesis of Finn Warlimont
%%=========================================================================

% Dimensions

length_total = 0.400;     %Barbell length [m]
length_middle = 0.200;    %Middle part length [m]
length_end = (length_total-length_middle)/2;


thickness_end = 0.2;    %Thickness of both ends [m]
thickness_middle = 0.1;   %Thickness of thin middle part [m]

% Material Constants

rho = 7850;      %Density [kg/m^3]
E = 210e9;       %Young's Modulus [N/m^2]
nu = 0.3;        %Contraction [-]

% Gravity

g = 9.81;       %Gravitational constant [m/s^2]

% Mesh-Resolution

NumberOfElementsX = 100;
NumberOfElementsY = 20;



%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY, E, nu, rho);

alpha = 0.00001;
beta = 0.0001;


D = alpha*K+beta*M;
%% Solve static FEM with boundary conditions

U_static = StaticFEM(K,ForceBoundaryCondition(NodeGrid,0),NodeGrid);
U_mass = StaticFEM(K,g*M*ones(size(U_static)),NodeGrid);
PlotDisplacement(U_static,NodeGrid,NodePosition)

%% Transient response

[t, U_dyn] = DynamicFEM(K,M,D,NodeGrid, [U_static/norm(U_static), U_mass/norm(U_mass)]);



