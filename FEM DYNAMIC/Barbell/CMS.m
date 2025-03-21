
function [M_tilde_ii,D_tilde_ii,K_tilde_ii,M_tilde_ie,D_tilde_ie,K_tilde_ie,V_cms] = CMS(K,M,D,NumberOfModes, AdditionalModes, NodeGrid)
    nno = NodeGrid(end,end);
    [~,boundNodes] = PositionBoundaryCondition(NodeGrid,0);

    % Reconfiguration of matrices
    % [ A_ee A_ei;
    %   A_ie A_ii]
    [K, K_ee, K_ei, K_ie, K_ii]= MatrixReconfiguration(K, boundNodes);
    [M, M_ee, M_ei, M_ie, M_ii]= MatrixReconfiguration(M, boundNodes);
    [D, D_ee, D_ei, D_ie, D_ii]= MatrixReconfiguration(D, boundNodes);


    % Eigenmodes of the free system
    if NumberOfModes == 0
        V_cms = [eye(length(boundNodes)); -K_ii\K_ie];
    else
        [Phi, Omega] = eigs(K_ii,M_ii,NumberOfModes,'smallestabs');
        Phi = [Phi, AdditionalModes];
        
        V_cms = [eye(length(boundNodes)) zeros(length(boundNodes),length(Phi(1,:)));
                 -K_ii\K_ie Phi];
    end
    
    % Calculate Matrices
    M_tilde = V_cms.'*M*V_cms;
    D_tilde = V_cms.'*D*V_cms;
    K_tilde = V_cms.'*K*V_cms;

    [K_tilde, K_tilde_ee, K_tilde_ei, K_tilde_ie, K_tilde_ii]= MatrixReconfiguration(K_tilde, 1:length(boundNodes));
    [M_tilde, M_tilde_ee, M_tilde_ei, M_tilde_ie, M_tilde_ii]= MatrixReconfiguration(M_tilde, 1:length(boundNodes));
    [D_tilde, D_tilde_ee, D_tilde_ei, D_tilde_ie, D_tilde_ii]= MatrixReconfiguration(D_tilde, 1:length(boundNodes));

end