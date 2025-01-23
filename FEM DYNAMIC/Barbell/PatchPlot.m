%% Visualisation

% defines vertices
tmp = NodePosition(:).';
U = U_dyn_modalFull(:,end);
verts = [tmp(1:2:end).'+U(1:2:end), tmp(2:2:end).'+U(2:2:end)];

% define faces
faces = [1, 2, 7, 6; ...
         2, 3, 8, 7; ...
         3, 4, 9, 8; ...
         4, 5, 10, 9; ...
         6, 7, 12, 11; ...
         7, 8, 13, 12; ...
         8, 9, 14, 13; ...
         9, 10, 15, 14];
verts = [tmp(1:2:end).'+U(1:2:end), tmp(2:2:end).'+ U(2:2:end)];
% create patch object
pObj = patch('vertices', verts, 'faces', faces, 'facecolor', 'none')


% for ii=1:1000:length(t_direct)
%     %  determine displacement    
%     U = U_dyn_direct(:,ii);
%     verts = [tmp(1:2:end).'+1000*U(1:2:end), tmp(2:2:end).'+ 1000*U(2:2:end)];
% 
% 
%     % update vertice position
%     set(pObj, 'vertices', verts);
% 
%     % update patch object
%     drawnow;
% 
% end


% v = u(1, :); % Extract the first frame of u
% v = reshape(v, 2, 15).'; % Reshape to 15 vertices of 2D coordinates
% f = [1, 2, 7, 6]; % Faces for the patch object
% f = [f; f+1; f+2; f+3]; % Extend faces
% f = [f; f+5]; % Further extend for visualization
% clf; % Clear current figure
% % Create the patch object for u
% pObj = patch('vertices', v, 'faces', f, 'FaceColor', 'cyan', 'EdgeColor', 'black');
% % Loop through each frame of u
% for ii = 1:2000:size(u, 1)
%     % Update the patch vertices
%     set(pObj, 'Vertices', reshape(u(ii, :), 2, 15).');
% 
%     % Clear previous h1 plot
%     hold on; % Retain current plot
%     if exist('h1_plot', 'var') && isvalid(h1_plot)
%         delete(h1_plot); % Remove previous h1 visualization
%     end
% 
%     % Extract h1 as a 2D vector from the current frame of u
%     h1 = [u(ii, 11), u(ii, 12)];
% 
%     % Plot h1 from origin to its endpoint
%     h1_plot = plot([0, h1(1)], [0, h1(2)], 'r-', 'LineWidth', 2); % Red line for h1
% 
%     % Clear previous timer text
%     if exist('timer_text', 'var') && isvalid(timer_text)
%         delete(timer_text); % Remove previous timer text
%     end
% 
%     % Add timer display
%     timer_text = text(-1, 3, sprintf('Time: %.2f s', t(ii)), 'FontSize', 12, 'Color', 'blue');
% 
%     % Set axis properties
%     axis equal;
%     xlim([-L0-1, L0+L+1]); % Adjust limits as necessary for your data
%     axis equal;    
%     grid on;
% 
%     % Update the plot
%     drawnow;
% end