%Gaussian Quadrature to generate stiffness and mass matrices

function [K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY, E, nu, rho)

TotalNumberOfElements = NumberOfElementsX*NumberOfElementsY;

    t = 1; 
    
    C = E/(1-nu^2)*[1 nu 0;...
                     nu 1 0;...
                     0 0 (1-nu)/2];

%    elementLoadVector = zeros(8,TotalNumberOfElements);
    KmatrixElement = zeros(8,8,TotalNumberOfElements);
    elementMassMatrix = zeros(8,8,TotalNumberOfElements);

    n = 2;      %Degree of quadrature
    [xi_vector, xi_weights] = GaussianQuadrature1D(2);
    [eta_vector, eta_weights] = GaussianQuadrature1D(2);


for e=1:TotalNumberOfElements

    dxdxi = 0.25*(-NodePositionTable(1)+NodePositionTable(3)+NodePositionTable(5)-NodePositionTable(7));
    dxdeta = 0.25*(-NodePositionTable(1)-NodePositionTable(3)+NodePositionTable(5)+NodePositionTable(7));
    dydxi = 0.25*(-NodePositionTable(2)+NodePositionTable(4)+NodePositionTable(6)-NodePositionTable(8));
    dydeta = 0.25*(-NodePositionTable(2)-NodePositionTable(4)+NodePositionTable(6)+NodePositionTable(8));
            
    invJ = [dxdxi dxdeta; dydxi dydeta];
    
    detJ = 1/det(invJ);

    for j = 1:n
        for i = 1:n
        
            B=B_matrix(xi_vector(i),eta_vector(j),J(1,1),J(1,2),J(2,1),J(2,2));
            KmatrixElement(:,:,e) = KmatrixElement(:,:,e)+B.'*C*B*t*detJ*xi_weights(i)*eta_weights(j);
            
            [N_temp,~,~] = ShapeFunctions(xi_vector(i),eta_vector(j));
            N = [N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4) 0;
                 0 N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4)];

            elementMassMatrix(:,:,e) = elementMassMatrix(:,:,e) + rho * (N.')*N*detJ*t*xi_weights(i)*eta_weights(j);
        end
    end
end
% Assembly
K = zeros(2*NumberOfNodes);
M = zeros(2*NumberOfNodes);

for e = 1:TotalNumberOfElements

    for i = 1:8
        for j =1:8
            K(NodeTable(e,i),NodeTable(e,j))=K(NodeTable(e,i),NodeTable(e,j))+KmatrixElement(i,j,e);
            M(NodeTable(e,i),NodeTable(e,j))=M(NodeTable(e,i),NodeTable(e,j))+elementMass(i,j,e); 
        end
    end
end
%% 


end
