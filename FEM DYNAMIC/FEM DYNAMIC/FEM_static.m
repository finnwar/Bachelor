


%Simulation eines einfachen Stabes mithilfe einer 2D-FEM

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
NumberOfElementsX = 15; % Unterteilung der X-Richtung
NumberOfElementsY = 5; % Unterteilung der Y-Richtung
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

%Displacement due to Force
U =K\f;
[n,~] = size(U);
U_x = zeros(n/2,1);
U_y = zeros(n/2,1);
% Displacement mode due to gravity
U_Gravity = K\(M*g*ones(size(f)));

%Eigenmdoes
[eigenVectors,eigenFrequencies] = eig(K,M);

%Damping Matrix via Rayleigh-Coefficients
alpha = 5;
beta = 0.005;

D = alpha * M + beta * K;

% Model Reduction
Phi = [U,U_Gravity,eigenVectors(:,1:4)];

M_tilde = Phi'*M*Phi;
K_tilde = Phi'*M*Phi;
D_tilde = Phi'*D*Phi;
f_tilde = Phi'*f;

% Transient Analysis


M_tilde_inv = inv(M_tilde);
%opt = odeset('maxstep',1e-9);
tspan = [0 1];

U_0 = U;
q_0 = Phi\U_0;

qdot_0 = zeros(length(q_0),1);
qddot_0 = zeros(length(q_0),1);

f_tilde_time = zeros(size(f_tilde));
Y0 = [q_0;qdot_0];
A = [zeros(size(K_tilde)) eye(size(K_tilde));
     -M_tilde_inv*K_tilde -M_tilde_inv*D_tilde];
b = [zeros(size(f_tilde_time));M_tilde_inv*f_tilde_time];

[t,Ysol] = ode45(@(t,Y) timeStepIntegration(t,Y,A,b), tspan, Y0);


U_sol = Phi * Ysol(:,1:length(Ysol(1,:))/2)';
figure
plot(t,U_sol(end/4,:))