function [StressFieldNodeGrid, StressFieldAbscissaGrid, Xgrid, Ygrid] = StressField(ElementStress,NodePositionTable,NodePositionGrid,nelX,nelY,xiVector,etaVector)
    
    StressFieldAbscissaGrid = zeros(2*nelX,2*nelY);
    Xgrid = zeros(2*nelX,2*nelY);
    Ygrid = zeros(2*nelX,2*nelY);

    for eX = 1:nelX
        for eY = 1:nelY
            e = eX + (eY-1)*nelX;
            StressFieldAbscissaGrid([2*eX-1 2*eX],[2*eY-1 2*eY]) = [ElementStress(e,1) ElementStress(e,2); ...
                                                        ElementStress(e,4) ElementStress(e,3)];
            N1 = ShapeFunctions(xiVector(1),etaVector(1));
            N2 = ShapeFunctions(xiVector(2),etaVector(1));
            N3 = ShapeFunctions(xiVector(2),etaVector(2));
            N4 = ShapeFunctions(xiVector(1),etaVector(2));
            
            AbscissaPos = [N1(1) 0 N1(2) 0 N1(3) 0 N1(4) 0; 
                           0 N1(1) 0 N1(2) 0 N1(3) 0 N1(4);
                           N2(1) 0 N2(2) 0 N2(3) 0 N2(4) 0; 
                           0 N2(1) 0 N2(2) 0 N2(3) 0 N2(4);
                           N3(1) 0 N3(2) 0 N3(3) 0 N3(4) 0; 
                           0 N3(1) 0 N3(2) 0 N3(3) 0 N3(4);
                           N4(1) 0 N4(2) 0 N4(3) 0 N4(4) 0; 
                           0 N4(1) 0 N4(2) 0 N4(3) 0 N4(4);];
            AbscissaPos = AbscissaPos * NodePositionTable(e,:).';
            Xgrid([2*eX-1 2*eX],[2*eY-1 2*eY])       = [AbscissaPos(1) AbscissaPos(3);
                                                        AbscissaPos(7) AbscissaPos(5)];
            Ygrid([2*eX-1 2*eX],[2*eY-1 2*eY])       = [AbscissaPos(2) AbscissaPos(3);
                                                        AbscissaPos(8) AbscissaPos(6)];   
        end
    end
    StressFieldInterpolant = scatteredInterpolant((Xgrid(:)),(Ygrid(:)),(StressFieldAbscissaGrid(:)));
    Xq = NodePositionGrid(1:2:end,:);
    Yq = NodePositionGrid(2:2:end,:);
    
    StressFieldNodeGrid = StressFieldInterpolant(Xq(:),Yq(:));
    StressFieldNodeGrid = reshape(StressFieldNodeGrid,size(Xq));
end