%% Visualisation
function PatchPlot(Title,U,t,Phi_vM,PhiX,PhiY,PhiXY,NodeGrid,NodePosition,NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle,isCenter)




% defines vertices
tmp = NodePosition(:).';

verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+U(2:2:end,1)];





% define faces
[~, ~, ~, ~, faces] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle);
% Calculate Stress at points defined in opt
if isCenter ==1
    MeanNodeStress = zeros(NumberOfElementsY*NumberOfElementsX);
else    
    MeanNodeStress = zeros(NodeGrid(end,end)/2,length(t));
end
for T = 1:length(t)
    MeanNodeStress(:,T) = sqrt((PhiX*U(:,T)).^2+(PhiY*U(:,T)).^2+abs((PhiX*U(:,T)).*(PhiY*U(:,T)))+3*(PhiXY*U(:,T)).^2);
end
verts = [tmp(1:2:end).'+U(1:2:end,1), tmp(2:2:end).'+ U(2:2:end,1)];
vertStress = MeanNodeStress(:,:,1);


% create patch object
figure
if isCenter ==1
    pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress(:,1),'FaceColor','flat')
else
    pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress(:,1),'FaceColor','interp')
end
daspect([1 1 1])
title(Title)
xlabel('x[m]')
ylabel('y[m]')
c = colorbar;
c.Label.String = 'Stress [Pa]';

axis('manual');
% set(gca,'ColorScale','log')
clim([0 max(max(vertStress))])
delta = 1;
if length(t) > 1000
    delta = round(length(t)/100);
end
axis([-0.25*length_end 1.2*(2*length_end+length_middle) -thickness_end thickness_end])
    t_delta = t(1:delta:length(t));
    iii=2;
    for ii=2:delta:length(t)
        pause(t_delta(iii)-t_delta(iii-1))
        %  determine displacement    
        verts = [tmp(1:2:end).'+1000*U(1:2:end,ii), tmp(2:2:end).'+ 1000*U(2:2:end,ii)];
        % update vertice position
        if isCenter==1
            set(pObj, 'vertices', verts, 'FaceVertexCData',MeanNodeStress(:,ii),'FaceColor','flat');
        else
            set(pObj, 'vertices', verts, 'FaceVertexCData',MeanNodeStress(:,ii));
        end
        % update patch object
        drawnow;
        iii=iii+1;
    end
    toc;
end
