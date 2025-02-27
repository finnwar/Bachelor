function [Phi_vonMises,Phi_X,Phi_Y,Phi_XY] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable, nu, E)
    % Stress modes evaluated at abscissa and extrapolated to node
    Phi_vonMises = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    
    Phi_X = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    Phi_Y = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    Phi_XY = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    % Stress modes evaluated directly at the noodes
    Phi_Xnode = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    Phi_Ynode = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    Phi_XYnode = zeros(NodeGrid(end,end)/2,NodeGrid(end,end));
    % Stress modes evaluated at the center
    % Phi_Xcenter = zeros(NodeTable(:,1))
    for i = 1:NodeGrid(end,end)
        unitLoad = zeros(NodeGrid(end,end),1);
        unitLoad(i) = 1;
    
        [vonMises,~, sigmaX,sigmaY,tauXY,sigmaXnode,sigmaYnode,tauXYnode,~,~,~] = StressCalculation(unitLoad,1,nu,E,NodePositionTable,NodeTable);
        Phi_vonMises(:,i) = StressField(vonMises,NodeTable,NodeGrid);
        Phi_X(:,i) = StressField(sigmaX,NodeTable,NodeGrid);
        Phi_Y(:,i) = StressField(sigmaY,NodeTable,NodeGrid);
        Phi_XY(:,i) = StressField(tauXY,NodeTable,NodeGrid);

        Phi_Xnode(:,i) = StressField(sigmaXnode,NodeTable,NodeGrid);
        Phi_Ynode(:,i) = StressField(sigmaYnode,NodeTable,NodeGrid);
        Phi_XYnode(:,i) = StressField(tauXYnode,NodeTable,NodeGrid);

    end
end