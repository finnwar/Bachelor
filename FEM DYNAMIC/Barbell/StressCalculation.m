function [vonMisesStressNode,vonMisesStressAbscissa, sigmaXnode,sigmaYnode,tauXYnode,sigmaXnodeDir,sigmaYnodeDir,tauXYnodeDir,...
    sigmaXCenter,sigmaYCenter,tauXYCenter]= StressCalculation(U,t,nu,E,NodePositionTable,NodeTable)

% Initialize arrays containing stress at abscissa
sigmaX = zeros(length(NodeTable(:,1)), 4, length(t));
sigmaY = zeros(length(NodeTable(:,1)), 4, length(t));
tauXY = zeros(length(NodeTable(:,1)), 4, length(t));
% Initialize arrays containing stress at node, extrapolated from abscissa
sigmaXnode = zeros(length(NodeTable(:,1)), 4, length(t));
sigmaYnode = zeros(length(NodeTable(:,1)), 4, length(t));
tauXYnode = zeros(length(NodeTable(:,1)), 4, length(t));
% Initialize arrays containing stress at node, evaluated at node
sigmaXnodeDir = zeros(length(NodeTable(:,1)), 4, length(t));
sigmaYnodeDir = zeros(length(NodeTable(:,1)), 4, length(t));
tauXYnodeDir = zeros(length(NodeTable(:,1)), 4, length(t));
% Initialize arrays containing stress at center
sigmaXCenter = zeros(length(NodeTable(:,1)), length(t));
sigmaYCenter = zeros(length(NodeTable(:,1)), length(t));
tauXYCenter = zeros(length(NodeTable(:,1)), length(t));
% Gauss abscissae
xi=GaussianQuadrature1D(2);
eta = xi;

xi_node = [-1 1];
eta_node = xi_node;

C = E/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];

transformationMatrix = [1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
                        -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
                        1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
                        -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)];
nodeHelp = [1 2 4 3];
for e = 1:length(NodeTable)
    for T = 1:length(t)        
        for i = 1:2
            for j= 1:2
                % Stress calculated at the Gauß-point and extrapolated to
                % node
                [~,dNdxi,dNdeta] = ShapeFunctions(xi(i),eta(j));
                J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
                detJ = det(J);
                invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
                B=B_matrix(xi(i), eta(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
                sigmaTemp = C*B*U(NodeTable(e,:),T);
                sigmaX(e,nodeHelp(i+(j-1)*2),T) = sigmaTemp(1);
                sigmaY(e,nodeHelp(i+(j-1)*2),T) = sigmaTemp(2);
                tauXY(e,nodeHelp(i+(j-1)*2),T)  = sigmaTemp(3);

                [~,dNdxi,dNdeta] = ShapeFunctions(xi_node(i),eta_node(j));
                J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
                detJ = det(J);
                invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
                B=B_matrix(xi_node(i), eta_node(j), invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
                sigmaTemp = C*B*U(NodeTable(e,:),T);
                sigmaXnodeDir(e,nodeHelp(i+(j-1)*2),T) = sigmaTemp(1);
                sigmaYnodeDir(e,nodeHelp(i+(j-1)*2),T) = sigmaTemp(2);
                tauXYnodeDir(e,nodeHelp(i+(j-1)*2),T)  = sigmaTemp(3);
            end
        end
        [~,dNdxi,dNdeta] = ShapeFunctions(0,0);
        J = [dNdxi;dNdeta]*[NodePositionTable(e,1:2:7).' NodePositionTable(e,2:2:8).'];
        detJ = det(J);
        invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
        B=B_matrix(0, 0, invJ(1,1), invJ(2,1), invJ(1,2), invJ(2,2));
        sigmaTemp = C*B*U(NodeTable(e,:),T);
        sigmaXCenter(e,T) = sigmaTemp(1);
        sigmaYCenter(e,T) = sigmaTemp(2);
        tauXYCenter(e,T)  = sigmaTemp(3);

        sigmaXnode(e,:,T)=(transformationMatrix*sigmaX(e,:,T).').';
        sigmaYnode(e,:,T)=(transformationMatrix*sigmaY(e,:,T).').';
        tauXYnode(e,:,T)=(transformationMatrix*tauXY(e,:,T).').';
    end
end

vonMisesStressAbscissa = sqrt(sigmaX.^2+sigmaY.^2+abs(sigmaX.*sigmaY)+3*(tauXY.^2));

vonMisesStressNode = sqrt(sigmaXnode.^2+sigmaYnode.^2+abs(sigmaXnode.*sigmaYnode)+3*(tauXYnode.^2));
end