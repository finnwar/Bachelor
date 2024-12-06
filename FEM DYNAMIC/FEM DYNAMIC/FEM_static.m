


%Simulation eines einfachen Stabes mithilfe einer 2D-FEM

%Materialeigenschaften

nu = 0.3;       % Querkontraktion
E = 210e9;      % Elastizitätsmodul [N/m^2]
rho = 7850;     % Dichte [kg/m^3]


C = E/(1-nu^2)*[1 nu 0;...
                nu 1 0;...
                0 0 (1-nu)/2];

%Maße des Stabes

L = 1;     %Länge[m]
H = 0.1;      %Höhe[m]

%Load Function

forceDistribution = @(u,v) ([0;0]);

%Auflösung
NumberOfElementsX = 10; % Unterteilung der X-Richtung
NumberOfElementsY = 3; % Unterteilung der Y-Richtung
TotalNumberOfElements = NumberOfElementsX*NumberOfElementsY;

%% 

NumberOfNodesX = NumberOfElementsX+1;
NumberOfNodesY = NumberOfElementsY+1;


NumberOfNodes = (NumberOfElementsX+1)*(NumberOfElementsY+1);


NodeGrid = zeros(2*NumberOfNodesX,NumberOfNodesY);

%Numbering of global nodes
for i=1:(2*NumberOfNodes)
    NodeGrid(ind2sub(size(NodeGrid.'),i))=i;
end
%  Global nodes of elements
%   Element progression:
%  x _____________
% y|  1,2   3,4   5,6   7,8
%  |  9,10  11,12 13,14 15,16
%  |  17,18 19,20 21,22 23,24
NodeTable = zeros(NumberOfElementsX*NumberOfElementsY,8);

for i = 1:(NumberOfElementsX)
    for j = 1:(NumberOfElementsY)
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),1) = NodeGrid(2*i-1,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),2) = NodeGrid(2*i,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),3) = NodeGrid(2*i+1,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),4) = NodeGrid(2*i+2,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),5) = NodeGrid(2*i+1,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),6) = NodeGrid(2*i+2,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),7) = NodeGrid(2*i-1,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),8) = NodeGrid(2*i,j+1);
        
    end
end

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
externalForce = [-1000000 NodeTable(end,end)];

%% 


dxidx = 2 / (L/NumberOfElementsX);
detadx = 0;
dxidy = 0;
detady = 2/(H/NumberOfElementsY);

dxdxi = 1/dxidx;
dxdeta = 0;
dydxi = 0;
dydeta = 1/detady;
% Jacobian des Elements (hier konstant wegen gewählter Geometrie
J = [dxdxi dydxi; dxdeta dydeta];
detJ =  det(J);           
% Element Steiffigkeitsmatrix

t = 1;

elementLoadVector= zeros(8,TotalNumberOfElements);
KmatrixElement=zeros(8,8,TotalNumberOfElements);
elementMass = zeros(8,8,TotalNumberOfElements);

for e = 1:TotalNumberOfElements

    % Gauss Quadratur
    n = 2; % Grad der Quadratur
    [xi_vector, xi_weights] = GaussianQuadrature1D(n);
    [eta_vector, eta_weights] = GaussianQuadrature1D(n);
    for j=1:n
        for i =1:n
            
            
            J = [dxdxi dydxi; dxdeta dydeta];
            
            B = B_matrix(xi_vector(i),eta_vector(j),dxidx,dxidy,detadx,detady);
            KmatrixElement(:,:,e) = KmatrixElement(:,:,e) + B.'*C*B*t*detJ*xi_weights(i)*eta_weights(j);

            [N_temp,~,~] = ShapeFunctions(xi_vector(i),eta_vector(j));
            N = [N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4) 0;
                 0 N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4)];
            elementLoadVector(:,e) = elementLoadVector(:,e)+N.'*forceDistribution(x(ind2sub(size(x),e)),y(ind2sub(size(y),e)))*t*detJ*xi_weights(i)*eta_weights(j);
        
            elementMass(:,:,e) = elementMass(:,:,e) + rho * (N.')*N*detJ*t*xi_weights(i)*eta_weights(j);
        end
    end
end
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


U =K\f;
[n,~] = size(U);
U_x = zeros(n/2,1);
U_y = zeros(n/2,1);
%% 


M_inv = inv(M);
%opt = odeset('maxstep',1e-9);
tspan = [0 1e-2];

U_0 = U;
Udot_0 = zeros(length(U_0),1);
Uddot_0 = zeros(length(U_0),1);

fTime = zeros(size(f));
Y0 = [U_0;Udot_0];
A = [zeros(size(K)) eye(size(K));
     -M*K zeros(size(K))];
b = [zeros(size(f));M_inv*fTime];

[t,Ysol] = ode45(@(t,Y) timeStepIntegration(t,Y,A,b), tspan, Y0);

figure
plot(t,Ysol(:,end/4))

%% 

