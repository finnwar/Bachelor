%Gaussian Quadrature to generate stiffness and mass matrices

function [K,M] = GaussianQuadrature(NodeTable, NodePositionTable, NumberOfElementsX, NumberOfElementsY,thickness, E, nu, rho)

TotalNumberOfElements = NumberOfElementsX*NumberOfElementsY;
NumberOfNodes = (NumberOfElementsX+1)*(NumberOfElementsY+1);
    thickness = 0.1; 
    
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

    


    for j = 1:n
        for i = 1:n
            
            %Evaluate Shapefunctions and their derivatives
            [N_temp,dNdxi,dNdeta] = ShapeFunctions(xi_vector(i),eta_vector(j));

            J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
            detJ = det(J);
            invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
            B=B_matrix(xi_vector(i), eta_vector(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
            KmatrixElement(:,:,e) = KmatrixElement(:,:,e)+B.'*C*B*thickness*detJ*xi_weights(i)*eta_weights(j);
            
            N = [N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4) 0;
                 0 N_temp(1) 0 N_temp(2) 0 N_temp(3) 0 N_temp(4)];

            elementMassMatrix(:,:,e) = elementMassMatrix(:,:,e) + rho * (N.')*N*detJ*thickness*xi_weights(i)*eta_weights(j);
        end
    end
end
% Assembly
K = zeros(2*NumberOfNodes);
M = zeros(2*NumberOfNodes);

for e = 1:TotalNumberOfElements

    for i = 1:8
        for j =1:8
            K(NodeTable(e,i), NodeTable(e,j))=K(NodeTable(e,i), NodeTable(e,j)) + KmatrixElement(i,j,e);
            M(NodeTable(e,i), NodeTable(e,j))=M(NodeTable(e,i), NodeTable(e,j)) + elementMassMatrix(i,j,e); 
        end
    end
end
%% 


end
