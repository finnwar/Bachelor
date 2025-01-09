%%Function to calculate the position of a given node

function [Pos, PosTable] = NodePositionCalc(NumberOfElementsX, NumberOfElementsY, length_middle,length_end, thickness_middle, thickness_end )
    

    Pos = zeros(2*(NumberOfElementsX+1),NumberOfElementsY+1);
    
    factors = GeometricParameterization(length_end, thickness_end, thickness_middle);


    for i = 1:(NumberOfElementsX+1)
        for j = 1:(NumberOfElementsY+1)
            X = ((2*length_end+length_middle)/NumberOfElementsX)*(i-1);
            
            Edge = EdgeFunction(factors, X, length_end, length_middle, thickness_middle);
            
            f = Edge(1);
            g = Edge(2);

            Y = (j-1)*(f-g)/NumberOfElementsY+g;
            Pos(2*i-1,j) = X;
            Pos(2*i, j) = Y;
        end
    end
    
    PosTable = zeros(NumberOfElementsX*NumberOfElementsY,8);

    for i = 1:(NumberOfElementsX)
    for j = 1:(NumberOfElementsY)
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),1) = Pos(2*i-1,j);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),2) = Pos(2*i,j);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),3) = Pos(2*i+1,j);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),4) = Pos(2*i+2,j);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),5) = Pos(2*i+1,j+1);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),6) = Pos(2*i+2,j+1);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),7) = Pos(2*i-1,j+1);
        PosTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),8) = Pos(2*i,j+1);
        
    end
end

end