function [sigma_node, sigma_abscissae] = StressCalculation(U,t,nu,E,NodePositionTable,NodeTable)
    elementDisplacement = zeros([size(NodeTable) length(t)]);
    for e = 1:length(NodeTable(:,1))
        for i = 1:length(NodeTable(1,:))
            elementDisplacement(e,i,:)= U(NodeTable(e,i),:);
        end
    end
    C = E/(1-nu^2)*[1 nu 0;...
        nu 1 0;...
        0 0 (1-nu)/2];
    sigma_node=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_abscissae=zeros(length(elementDisplacement(:,1)),4,length(t));
    
    [xi_vector,  ~] = GaussianQuadrature1D(2);
    [eta_vector, ~] = GaussianQuadrature1D(2);
    AbscissaPositionTable = zeros(size(NodePositionTable));
    xi_nodes  = [-1 1];
    eta_nodes = [-1 1];
    for T = 1:length(t)
        for e=1:length(elementDisplacement(:,1))
            for j = 1:2
                for i = 1:2
                    [~,dNdxi,dNdeta] = ShapeFunctions(xi_nodes(i),eta_nodes(j));
                    J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
                    detJ = det(J);
                    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
    
                    B = B_matrix(xi_nodes(i),eta_nodes(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
    
                    sigma_node(e, i+2*(j-1),T) = vecnorm(C*B*elementDisplacement(e,:,T).');

                    [~,dNdxi,dNdeta] = ShapeFunctions(xi_vector(i),eta_vector(j));
                    J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
                    detJ = det(J);
                    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
                    B=B_matrix(xi_vector(i), eta_vector(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));

                    sigma_abscissae(e, i+2*(j-1),T) = vecnorm(C*B*elementDisplacement(e,:,T).');
                end
            end
        end
    end

end