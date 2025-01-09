%% Function to generate the complete Mesh Information

function [NodeGrid, NodeTable, NodePosition, NodePositionTable] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)

    [NodeGrid, NodeTable] = NodeGridGenerator(NumberOfElementsX, NumberOfElementsY);

    [NodePosition, NodePositionTable] = NodePositionCalc(NumberOfElementsX, NumberOfElementsY, ...
                    length_middle, length_end, thickness_middle, thickness_end);

    



end