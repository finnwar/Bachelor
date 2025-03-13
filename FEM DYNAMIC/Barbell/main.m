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
maxU = zeros(20,1);
avgU=maxU;
for i =1:20
NumberOfElementsX = 2*i;
NumberOfElementsY = i;

%% Mesh generation
[NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX, NumberOfElementsY, length_end, length_middle, thickness_end, thickness_middle);

%% Gaussian Quadrature to receive K, and M Matrices

[K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho);

alpha = 0;
beta = 1e-3;

D = alpha*K+beta*M;


% [Phi_vM,PhiX,PhiY,PhiXY,PhiXnode,PhiYnode,PhiXYnode,PhiXcenter,PhiYcenter,PhiXYcenter] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable,nu,E);

%% Transient response

[t_dir, U_dyn_dir] = DynamicFEM(K,M,D,NodeGrid);

endU = sqrt(U_dyn_dir(end-1,:).^2+U_dyn_dir(end,:).^2);
maxU(i)=max(endU);
avgU(i) = mean(endU);
end
%%
deltaMaxU=zeros(19,1);
relDeltaMaxU= deltaMaxU;

deltaAvgU=zeros(19,1);
relDeltaAvgU= deltaMaxU;

for i = 1:length(maxU)-1
    deltaMaxU(i) = abs(maxU(i+1)-maxU(i));
    relDeltaMaxU(i)=maxU(i+1)/maxU(i);
    
    deltaAvgU(i)=abs(avgU(i+1)-avgU(i));
    relDeltaAvgU(i)=avgU(i+1)/avgU(i);
end