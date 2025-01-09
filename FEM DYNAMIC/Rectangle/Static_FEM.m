function [K, M, f, U] = Static_FEM(NumberOfElementsX, NumberOfElementsY)
%Materialeigenschaften

nu = 0.3;       % Querkontraktion
E = 210e9;      % Elastizitätsmodul [N/m^2]
rho = 7850;     % Dichte [kg/m^3]
g = -9.81;       % Gravitationskonstante [m/s^2]


C = E/(1-nu^2)*[1 nu 0;...
                nu 1 0;...
                0 0 (1-nu)/2];

%Maße des Stabes

L = 1;     %Länge[m]
H = 0.1;      %Höhe[m]

%Load Function

forceDistribution = @(u,v) ([0;0]);

%Auflösung

TotalNumberOfElements = NumberOfElementsX*NumberOfElementsY;

%% 

NumberOfNodesX = NumberOfElementsX+1;
NumberOfNodesY = NumberOfElementsY+1;


NumberOfNodes = (NumberOfElementsX+1)*(NumberOfElementsY+1);

[NodeGrid, NodeTable] = NodeGridGenerator(NumberOfNodesX, NumberOfNodesY, NumberOfNodes, NumberOfElementsX, NumberOfElementsY);

%Grid
[x,y] = ndgrid(linspace(0,L,NumberOfNodesX),...
               linspace(0,H,NumberOfNodesY));

%% 

% Boundary Conditions

% Dirichlet Boundary Condition
% for a fixed rod
% displacementBoundaryCondition(value, node)

% Initialization
displacementBoundaryCondition = zeros(2*NumberOfNodesY,2);
% node
for i = 1:NumberOfNodesY
    displacementBoundaryCondition(2*i-1,2)=NodeGrid(1,i);
    displacementBoundaryCondition(2*i,2)=NodeGrid(2,i);
end

% Neumann Boundary Condition
externalForce = [1000 2*NumberOfNodes];

[elementLoadVector, KmatrixElement, elementMass] = GaussQuadratur(L, NumberOfElementsX, H, NumberOfElementsY, TotalNumberOfElements, C, forceDistribution, x, y, rho);
%% 

% Assembly
K = zeros(2*NumberOfNodes);
M = zeros(2*NumberOfNodes);
f = zeros(2*NumberOfNodes,1);
for e = 1:TotalNumberOfElements
    f(NodeTable(e,:))=f(NodeTable(e,:))+elementLoadVector(:,e);
    for i = 1:8
        for j =1:8
            K(NodeTable(e,i),NodeTable(e,j))=K(NodeTable(e,i),NodeTable(e,j))+KmatrixElement(i,j,e);
            M(NodeTable(e,i),NodeTable(e,j))=M(NodeTable(e,i),NodeTable(e,j))+elementMass(i,j,e); 
        end
    end
end
%% 

% Application of boundary conditions
% Neumann

[n,~] = size(externalForce);

for i = 1:n
    f(externalForce(i,2))=f(externalForce(i,2))+externalForce(i,1);
end
% Dirichlet

K(:,displacementBoundaryCondition(:,2))=[];
K(displacementBoundaryCondition(:,2),:)=[];

M(:,displacementBoundaryCondition(:,2))=[];
M(displacementBoundaryCondition(:,2),:)=[];


f(displacementBoundaryCondition(:,2)) = [];
%% 

%Displacement due to Force
U =K\f;
[n,~] = size(U);
U_x = zeros(n/2,1);
U_y = zeros(n/2,1);
% Displacement mode due to gravity
U_Gravity = K\(M*g*ones(size(f)));
end