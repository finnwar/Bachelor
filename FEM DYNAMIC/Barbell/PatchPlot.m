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
axis('manual');
axis([-0.25*length_end 1.2*(2*length_end+length_middle) -thickness_end thickness_end])
for ii=1:length(t)
    %  determine displacement    
    
    verts = [tmp(1:2:end).'+1000*U(1:2:end,ii), tmp(2:2:end).'+ 1000*U(2:2:end,ii)];
    for y = 2:(NumberOfElementsY-1)
        for x = 2:(NumberOfElementsX-1)
            e = x+(y-1)*NumberOfElementsX;
            vertStress(NodeTable(e,2:2:8)/2) = [0.25*(nodeStress(e,1,ii)+nodeStress(e-1-NumberOfElementsX,3,ii)+nodeStress(e-NumberOfElementsX,4,ii))+nodeStress(e-1,2,ii);
                                                0.25*(nodeStress(e,2,ii)+nodeStress(e-NumberOfElementsX,3,ii)+nodeStress(e+1-NumberOfElementsX,4,ii))+nodeStress(e+1,1,ii);
                                                0.25*(nodeStress(e,3,ii)+nodeStress(e+1,4,ii)+nodeStress(e+1+NumberOfElementsX,1,ii))+nodeStress(e+NumberOfElementsX,2,ii);
                                                0.25*(nodeStress(e,4,ii)+nodeStress(e-1+NumberOfElementsX,2,ii)+nodeStress(e+NumberOfElementsX,1,ii))+nodeStress(e-1,3,ii);];
        end
    end
    vertStress(NodeTable(1:NumberOfElementsX,[2 4])/2) = [0.5*(nodeStress(2:NumberOfElementsX,1,ii)+nodeStress(1:(NumberOfElementsX-1),2,ii));
                                                          0.5*(nodeStress(1:NumberOfElementsX,1,ii)+nodeStress(2:(NumberOfElementsX),2,ii))]
    % update vertice position
    set(pObj, 'vertices', verts, 'FaceVertexCData',vertStress);

    % update patch object
    drawnow limitrate;
    if ii < length(t)
        pause(t(ii+1)-t(ii))
    end
end
end
