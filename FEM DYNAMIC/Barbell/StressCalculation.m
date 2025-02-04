function [sigma_node, sigma_abscissae,sigma_nodeX,sigma_nodeY,sigma_nodeXY,sigma_abscissaeX,sigma_abscissaeY,sigma_abscissaeXY]...
                                                = StressCalculation(U,t,nu,E,NodePositionTable,NodeTable)
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
    sigma_nodeX=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_nodeY=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_nodeXY=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_abscissae=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_abscissaeX=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_abscissaeY=zeros(length(elementDisplacement(:,1)),4,length(t));
    sigma_abscissaeXY=zeros(length(elementDisplacement(:,1)),4,length(t));
    [xi_vector,  ~] = GaussianQuadrature1D(2);
    [eta_vector, ~] = GaussianQuadrature1D(2);

    for T = 1:length(t)
        for e=1:length(elementDisplacement(:,1))
            for j = 1:2
                for i = 1:2

                    [~,dNdxi,dNdeta] = ShapeFunctions(xi_vector(i),eta_vector(j));
                    J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
                    detJ = det(J);
                    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
                    B=B_matrix(xi_vector(i), eta_vector(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
                    temp = (C*B*elementDisplacement(e,:,T).');
                    sigma_abscissaeX(e, i+2*(j-1),T) = temp(1);
                    sigma_abscissaeY(e, i+2*(j-1),T) = temp(2);
                    sigma_abscissaeXY(e, i+2*(j-1),T) = temp(3);
                end
            end
        end
        
        for e=1:length(elementDisplacement(:,1))
            tempX = ([1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
                                 -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
                                 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
                                 -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)]*sigma_abscissaeX(e,:,T).').';
             tempY = ([1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
                                 -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
                                 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
                                 -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)]*sigma_abscissaeY(e,:,T).').';
             tempXY = ([1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
                                 -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
                                 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
                                 -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)]*sigma_abscissaeXY(e,:,T).').';
            sigma_nodeX(e,:,T) = tempX.';
            sigma_nodeY(e,:,T) = tempY.';
            sigma_nodeXY(e,:,T) = tempXY.';
            sigma_node(e,:,T) = sqrt(tempX.^2+tempY.^2-tempX.*tempY+3*tempXY.^2);
        end
    end
    sigma_abscissae = sqrt(sigma_abscissaeX.^2+sigma_abscissaeY.^2-sigma_abscissaeX.*sigma_abscissaeY+3*sigma_abscissaeXY.^2);
end