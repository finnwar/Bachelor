
function [V_cms,Phi,M_tilde,D_tilde,K_tilde,invK_iiKie] = ModalReduction(K,M,D,NumberOfModes, AdditionalModes, NodeGrid)
    nno = NodeGrid(end,end);


    [~,boundNodes] = PositionBoundaryCondition(NodeGrid,0);
    nFreeNodes = nno-length(boundNodes);
    nBroundNodes = length(boundNodes);
    
    %     Reconfiguration of matrices
    % [ A_ee A_ei;
    %   A_ie A_ii]
    [K, K_ee, K_ei, K_ie, K_ii]= MatrixReconfiguration(K, boundNodes);
    [M, M_ee, M_ei, M_ie, M_ii]= MatrixReconfiguration(M, boundNodes);
    [D, D_ee, D_ei, D_ie, D_ii]= MatrixReconfiguration(D, boundNodes);
    
    % Eigenmodes of the free system

    [Phi, ~] = eigs(K_ii,M_ii,NumberOfModes,'smallestabs');
    % 
    % [Omega, ind] = sort(Omega);
    % Phi = Phi(:, ind);
    Phi = [Phi(:,1:NumberOfModes), AdditionalModes];

   
    % Calculate Matrices

    invK_iiKie = K_ii\K_ie;
    V_cms = [eye(nBroundNodes) zeros(nBroundNodes,length(Phi(1,:)));
             -invK_iiKie Phi];


    M_tilde = V_cms.'*M*V_cms;
    D_tilde = V_cms.'*D*V_cms;
    K_tilde = V_cms.'*K*V_cms;
end