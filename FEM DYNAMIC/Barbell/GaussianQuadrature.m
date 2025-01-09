%Gaussian Quadrature to generate stiffness and mass matrices

function [K,M] = GaussianQuadrature(NodeGrid, NodeTable, NodePosition, NodePositionTable, NumberOfElementsX, NumberOfElementsY)

TotalNumberOfElements = NumberOfElementsX*NumberOfElementsY;

for e=1:TotalNumberOfElements
    n = 2;      %Degree of quadrature
    [xi_vector, xi_weights] = GaussianQuadrature1D(2);
    [eta_vector, eta_weights] = GaussianQuadrature1D(2);

    dxdxi = 0.25*(-NodePositionTable(1)+NodePositionTable(3)+NodePositionTable(5)-NodePositionTable(7));
    dxdeta = 0.25*(-NodePositionTable(1)-NodePositionTable(3)+NodePositionTable(5)+NodePositionTable(7));
    dydxi = 0.25*(-NodePositionTable(2)+NodePositionTable(4)+NodePositionTable(6)-NodePositionTable(8));
    dydeta = 0.25*(-NodePositionTable(2)-NodePositionTable(4)+NodePositionTable(6)+NodePositionTable(8));
            
    
    

    for j = 1:n
        for i = 1:n




end


