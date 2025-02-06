%% Visualisation
function PatchPlot(Title,U,t,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle)




% defines vertices
tmp = NodePosition(:).';

verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+U(2:2:end,1)];





% define faces
[~, ~, ~, ~, faces] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle);
% Calculate Stress at points defined in opt

MeanNodeStress = zeros(NodeGrid(end,end)/2,length(t));
% for T = 1:length(t)
%     MeanNodeStress(:,T) = Phi_vM*U(:,T);
% end
for T = 1:length(t)
    MeanNodeStress(:,T) = sqrt((PhiX*U(:,T)).^2+(PhiY*U(:,T)).^2+abs((PhiX*U(:,T)).*(PhiY*U(:,T)))+3*(PhiXY*U(:,T)).^2);
end
verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+ U(2:2:end,1)];
vertStress = MeanNodeStress(:,:,1);


% create patch object
figure
pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress(:,1),'FaceColor','interp')
title(Title)
xlabel('x')
ylabel('y')
c.Label.String = 'Stress [N/mm^2]';
colorbar
axis('manual');
clim([0 max(max(vertStress))])

axis([-0.25*length_end 1.2*(2*length_end+length_middle) -thickness_end thickness_end])
    for ii=2:length(t)
        pause(t(ii)-t(ii-1))
        %  determine displacement    
        verts = [tmp(1:2:end).'+1000*U(1:2:end,ii), tmp(2:2:end).'+ 1000*U(2:2:end,ii)];
        
        % update vertice position
        set(pObj, 'vertices', verts, 'FaceVertexCData',MeanNodeStress(:,ii));
        % update patch object
        drawnow limitrate;
    end
end
