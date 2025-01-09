function [NodeGrid, NodeTable] = NodeGridGenerator(NumberOfElementsX, NumberOfElementsY)

NumberOfNodesX = NumberOfElementsX+1;
NumberOfNodesY = NumberOfElementsY+1;
NumberOfNodes = NumberOfNodesX*NumberOfNodesY;

NodeGrid = zeros(2*NumberOfNodesX,NumberOfNodesY);

%Numbering of global nodes
for i=1:(2*NumberOfNodes)
    NodeGrid(ind2sub(size(NodeGrid.'),i))=i;
end
%  Global numbering of degrees of freedom
%  Progression:
%  x _____________
% y|  1,2   3,4   5,6   7,8
%  |  9,10  11,12 13,14 15,16
%  |  17,18 19,20 21,22 23,24
NodeTable = zeros(NumberOfElementsX*NumberOfElementsY,8);

for i = 1:(NumberOfElementsX)
    for j = 1:(NumberOfElementsY)
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),1) = NodeGrid(2*i-1,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),2) = NodeGrid(2*i,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),3) = NodeGrid(2*i+1,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),4) = NodeGrid(2*i+2,j);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),5) = NodeGrid(2*i+1,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),6) = NodeGrid(2*i+2,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),7) = NodeGrid(2*i-1,j+1);
        NodeTable(sub2ind([NumberOfElementsX,NumberOfElementsY],i,j),8) = NodeGrid(2*i,j+1);
        
    end
end
end