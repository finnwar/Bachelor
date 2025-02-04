%% Visualisation
function PatchPlot(Title,U,t,ElementStress,NodePosition,NodePositionTable,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle,opt)




% defines vertices
tmp = NodePosition(:).';

verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+U(2:2:end,1)];





% define faces
[~, ~, ~, ~, faces] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle);
% Calculate Stress at points defined in opt
if opt == 'abscissa'
    eta = GaussianQuadrature1D(2);
    xi = GaussianQuadrature1D(2);
elseif opt == 'node'
    eta = [-1 1];
    xi = [-1 1];
end

StressFieldNodeGrid = zeros([size(NodePosition(1:2:end,:)) length(t)]);
for T = 1:length(t)
    StressFieldNodeGrid(:,:,T) = StressField(ElementStress(:,:,T),NodePositionTable,NodePosition,NumberOfElementsX,NumberOfElementsY,xi,eta);
end
verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+ U(2:2:end,1)];
vertStress = StressFieldNodeGrid(:,:,1);


% create patch object
figure
pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress(:),'FaceColor','interp')
title(Title)
xlabel('x')
ylabel('y')
colorbar
axis('manual');
axis([-0.25*length_end 1.2*(2*length_end+length_middle) -thickness_end thickness_end])
for ii=1:length(t)
    %  determine displacement    
    
    verts = [tmp(1:2:end).'+1000*U(1:2:end,ii), tmp(2:2:end).'+ 1000*U(2:2:end,ii)];
    
    vertStress = StressFieldNodeGrid(:,:,ii);

    % update vertice position
    set(pObj, 'vertices', verts, 'FaceVertexCData',vertStress(:));

    % update patch object
    drawnow limitrate;
    if ii < length(t)
        pause(t(ii+1)-t(ii))
    end
end
end
