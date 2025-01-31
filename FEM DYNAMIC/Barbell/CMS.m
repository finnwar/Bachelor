
function [M_ii,D_ii,K_ii,K_ie,B_M,B_D,B_K,Phi,Omega] = CMS(K,M,D,NumberOfModes, AdditionalModes, NodeGrid)
    nno = NodeGrid(end,end);
    [~,boundNodes] = PositionBoundaryCondition(NodeGrid,0);

    % Reconfiguration of matrices
    % [ A_ee A_ei;
    %   A_ie A_ii]
    [K, K_ee, K_ei, K_ie, K_ii]= MatrixReconfiguration(K, boundNodes);
    [M, M_ee, M_ei, M_ie, M_ii]= MatrixReconfiguration(M, boundNodes);
    [D, D_ee, D_ei, D_ie, D_ii]= MatrixReconfiguration(D, boundNodes);
    K_ii_inv = inv(K_ii);

    % Eigenmodes of the free system

    [Phi, Omega] = eigs(K_ii,M_ii,NumberOfModes,'smallestabs');
    Phi = [Phi(:,1:NumberOfModes), AdditionalModes];

    % Calculate Matrices
    B_M = M_ie*eye(length(boundNodes))-M_ii*K_ii_inv*K_ie;
    B_D = D_ie*eye(length(boundNodes))-D_ii*K_ii_inv*K_ie;
    B_K = K_ie*eye(length(boundNodes))-K_ii*K_ii_inv*K_ie;
end