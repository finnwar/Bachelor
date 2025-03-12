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
NumberOfElementsY = 12;

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
try
tic;
[t_cms1, U_dyn_cms1] = DynamicCMSFEM(K,M,D,NodeGrid,1,[]);
toc;
catch
    disp("CMS 1 failed.")
end

%%
try
tic;
[t_cms5, U_dyn_cms5] = DynamicCMSFEM(K,M,D,NodeGrid,5,[]);
toc;
catch
    disp("CMS 5 failed.")
end
%%
try
tic;
[t_cms20, U_dyn_cms20] = DynamicCMSFEM(K,M,D,NodeGrid,20,[]);
toc;
catch
    disp("CMS 20 failed.")
end

%%
try

tic;
[t_cms30, U_dyn_cms30] = DynamicCMSFEM(K,M,D,NodeGrid,30,[]);
toc;
catch
    disp("CMS 30 failed.")
end
%%
try
tic;
[t_cms50, U_dyn_cms50] = DynamicCMSFEM(K,M,D,NodeGrid,50,[]);
toc;
catch
    disp("CMS 50 failed.")
end
%%
try
tic;
[t_cms80, U_dyn_cms80] = DynamicCMSFEM(K,M,D,NodeGrid,80,[]);
toc;
catch
    disp("CMS 80 failed.")
end
%%
try
tic;
[t_cms100, U_dyn_cms100] = DynamicCMSFEM(K,M,D,NodeGrid,100,[]);
toc;
catch
    disp("CMS 100 failed.")
end

%%
try
tic;
[t_cms120, U_dyn_cms120] = DynamicCMSFEM(K,M,D,NodeGrid,120,[]);
toc;
catch
    disp("CMS 120 failed.")
end
%%
try
tic;
[t_cms150, U_dyn_cms150] = DynamicCMSFEM(K,M,D,NodeGrid,150,[]);
toc;
catch
    disp("CMS 150 failed.")
end
%%
try
tic;
[t_cm200, U_dyn_cms200] = DynamicCMSFEM(K,M,D,NodeGrid,200,[]);
toc;
catch
    disp("CMS 200 failed.")
end
%%
try
tic;
[t_cms250, U_dyn_cms250] = DynamicCMSFEM(K,M,D,NodeGrid,250,[]);
toc;
catch
    disp("CMS 250 failed.")
end
%%
try
tic;
[t_cms270, U_dyn_cms270] = DynamicCMSFEM(K,M,D,NodeGrid,270,[]);
toc;
catch
    disp("CMS 270 failed.")
end
%%
try
tic;
[t_cms300, U_dyn_cms300] = DynamicCMSFEM(K,M,D,NodeGrid,300,[]);
toc;
catch
    disp("CMS 300 failed.")
end
%%
try
tic;
[t_cms400, U_dyn_cms400] = DynamicCMSFEM(K,M,D,NodeGrid,400,[]);
toc;
catch
    disp("CMS 400 failed.")
end
%%
try
tic;
[t_cms500, U_dyn_cms500] = DynamicCMSFEM(K,M,D,NodeGrid,500,[]);
toc;
catch
    disp("CMS 500 failed.")
end
%%
try
tic;
[t_cmsFull, U_dyn_cmsFull] = DynamicCMSFEM(K,M,D,NodeGrid,1000,[]);
toc;
catch
    disp("CMS Full failed.")
end

% %% Visualisation 
% %% Nodal Approach
% PatchPlot('Nodal Approach',U_dyn_dir,t_dir,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('CMS 1',U_dyn_cms1,t_cms1,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
%% 
% %% 
% 
% PatchPlot('CMS 20',U_dyn_cms20,t_cms20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('Modal 20',U_dyn_mod20,t_mod20,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% 
% %% 
% 
% PatchPlot('Modal 100',U_dyn_mod100,t_mod100,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('CMS 100',U_dyn_cms100,t_cms100,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('Modal 270',U_dyn_mod270,t_mod270,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('CMS 270',U_dyn_cms270,t_cms270,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)
% %% 
% 
% PatchPlot('Static',U_static,[1],Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
%                                                              length_end, length_middle, thickness_end, thickness_middle)








