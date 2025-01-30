%% Visualisation
function PatchPlot(Title,U,t,nodeStress,NodePosition,NodeTable,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)




% defines vertices
tmp = NodePosition(:).';

verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+U(2:2:end,1)];





% define faces
[~, ~, ~, ~, faces] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle);
% Calculate Stress at vertices

vertStress = zeros(length(verts(:,1)),1);
for e = 1:(NumberOfElementsY*NumberOfElementsX)
        vertStress(NodeTable(e,2:2:8)/2) = nodeStress(e,:,1).';
end
verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+ U(2:2:end,1)];



% create patch object
figure
pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress,'FaceColor','interp')
title(Title)
xlabel('x')
ylabel('y')
colorbar
for ii=1:length(t)
    %  determine displacement    
    
    verts = [tmp(1:2:end).'+1000*U(1:2:end,ii), tmp(2:2:end).'+ 1000*U(2:2:end,ii)];
    
    
    for e = 1:(NumberOfElementsY*NumberOfElementsX)
        vertStress(NodeTable(e,2:2:8)/2) = nodeStress(e,:,ii).';
    end

    % update vertice position
    set(pObj, 'vertices', verts, 'FaceVertexCData',vertStress);

    % update patch object
    drawnow limitrate;
    if ii < length(t)
        pause(t(ii+1)-t(ii))
    end
end
end
