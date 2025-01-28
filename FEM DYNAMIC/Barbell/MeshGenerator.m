%% Function to generate the complete Mesh Information

function [DoFGrid, DoFTable, NodePosition, NodePositionTable, vertices] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

    [DoFGrid, DoFTable] = NodeGridGenerator(NumberOfElementsX, NumberOfElementsY);

    [NodePosition, NodePositionTable] = NodePositionCalc(NumberOfElementsX, NumberOfElementsY, ...
                    length_middle, length_end, thickness_middle, thickness_end);

    vertexGrid = zeros(NumberOfElementsX+1,NumberOfElementsY+1);
    for i = 1:length(vertexGrid(:,1))
        for j = 1:length(vertexGrid(1,:))
            vertexGrid(i,j) = i+(j-1)*(NumberOfElementsX+1);
        end
    end

    vertices = zeros(NumberOfElementsX*NumberOfElementsY,4);

    for i = 1:(NumberOfElementsX)
        for j = 1:(NumberOfElementsY)
            vertices(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),1) = vertexGrid(i,j);
            vertices(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),2) = vertexGrid(i+1,j);
            vertices(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),3) = vertexGrid(i+1,j+1);
            vertices(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),4) = vertexGrid(i,j+1);

        end
    end


end