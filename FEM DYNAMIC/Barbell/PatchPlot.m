%% Visualisation

% defines vertices
tmp = NodePosition(:).';
U = U_static(:,end);
verts = [tmp(1:2:end).'+U(1:2:end), tmp(2:2:end).'+U(2:2:end)];





% define faces
[~, ~, ~, ~, faces] = MeshGenerator(NumberOfElementsX,NumberOfElementsY, ...
                                                             length_end, length_middle, thickness_end, thickness_middle);
% Calculate Stress at vertices

vertStress = zeros(length(verts(:,1)),1);
for e = 1:(NumberOfElementsY*NumberOfElementsX)
        vertStress(NodeTable(e,2:2:8)/2) = nodeStress(e,:,1).';
end
verts = [tmp(1:2:end).'+U(1:2:end), tmp(2:2:end).'+ U(2:2:end)];



% create patch object
pObj = patch('vertices', verts, 'faces', faces, 'FaceVertexCData',vertStress,'FaceColor','interp')



%% Animation the Direct Solver
for ii=1:length(t_dir)
    %  determine displacement    
    U = U_dyn_dir(:,ii);
    verts = [tmp(1:2:end).'+1000*U(1:2:end), tmp(2:2:end).'+ 1000*U(2:2:end)];
    
    
    for e = 1:(NumberOfElementsY*NumberOfElementsX)
        vertStress(NodeTable(e,2:2:8)/2) = nodeStress_dir(e,:,ii).';
    end

    % update vertice position
    set(pObj, 'vertices', verts, 'FaceVertexCData',vertStress);

    % update patch object
    drawnow;
    pause(t_dir(ii+1)-t_dir(ii))
end

%% Animation of the CMS Solver

for ii=1:length(t_mod)
    %  determine displacement    
    U = U_dyn_dir(:,ii);
    verts = [tmp(1:2:end).'+1000*U(1:2:end), tmp(2:2:end).'+ 1000*U(2:2:end)];
    
    
    for e = 1:(NumberOfElementsY*NumberOfElementsX)
        vertStress(NodeTable(e,2:2:8)/2) = nodeStress_mod(e,:,ii).';
    end

    % update vertice position
    set(pObj, 'vertices', verts, 'FaceVertexCData',vertStress);

    % update patch object
    drawnow;
    pause(t_dir(ii+1)-t_dir(ii))
end

