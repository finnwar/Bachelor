function [Phi_vonMises,Phi_X,Phi_Y,Phi_XY] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable, nu, E)
    Phi_vonMises = zeros(NodeGrid(end,end)/2);
    Phi_X = zeros(NodeGrid(end,end)/2);
    Phi_Y = zeros(NodeGrid(end,end)/2);
    Phi_XY = zeros(NodeGrid(end,end)/2);

    for i = 1:NodeGrid(end,end)
        unitLoad = zeros(NodeGrid(end,end),1);
        unitLoad(i) = 1;
    
        [vonMises,~, sigmaX,sigmaY,tauXY] = StressCalculation(unitLoad,1,nu,E,NodePositionTable,NodeTable);
        Phi_vonMises(:,i) = StressField(vonMises,NodeTable,NodeGrid);
        Phi_X(i,:) = StressField(sigmaX,NodeTable,NodeGrid);
        Phi_Y(i,:) = StressField(sigmaY,NodeTable,NodeGrid);
        Phi_XY(i,:) = StressField(tauXY,NodeTable,NodeGrid);
    end
end