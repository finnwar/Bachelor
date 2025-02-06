function [Phi_vonMises,Phi_X,Phi_Y,Phi_XY] = StressModeCalculation(NodeGrid,NodeTable,NodePositionTable, nu, E)
    Phi_vonMises = zeros(NodeGrid(end,end));
    Phi_X = zeros(NodeGrid(end,end));
    Phi_Y = zeros(NodeGrid(end,end));
    Phi_XY = zeros(NodeGrid(end,end));

    for i = 1:NodeGrid(end,end)
        unitLoad = zeros(NodeGrid(end,end),1);
        unitLoad(i) = 1;
    
        [vonMises,~, sigmaX,sigmaY,tauXY] = StressCalculation(uniLoad,1,nu,E,NodePositionTable,NodeTable);
        Phi_vonMises(:,i) = StressField(vonMises);
        Phi_X(i,:) = StressField(sigmaX);
        Phi_Y(i,:) = StressField(sigmaY);
        Phi_XY(i,:) = StressField(tauXY);
    end
end