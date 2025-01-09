function [elementLoadVector, KmatrixElement, elementMass] = GaussQuadratur(L, NumberOfElementsX, H, NumberOfElementsY, TotalNumberOfElements, C, forceDistribution, x, y, rho)
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
end