%% Calculates the average stress at the nodes
function [NodeStressMean] = StressField(ElementStress,NodeTable,NodeGrid)
    NodeStressMean = zeros(NodeGrid(end,end)/2,1);
    
    NodeCounter = zeros(NodeGrid(end,end)/2,1);       % Counts how many elements come together at a node
    for e = 1:length(NodeTable(:,1))
        NodeStressMean(NodeTable(e,2:2:end)/2) = NodeStressMean(NodeTable(e,2:2:end)/2) + ElementStress(e,:).';
        NodeCounter(NodeTable(e,2:2:end)/2) = NodeCounter(NodeTable(e,2:2:end)/2) + ones(4,1);
    end

    NodeStressMean = NodeStressMean./NodeCounter;
end