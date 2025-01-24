function [Omega, Phi, D_tilde, D_bar, M_bar, invK_FFK_FB] = ModalReduction(K,M,D,NumberOfModes, AdditionalModes, NodeGrid)
    
    [~,boundNodes] = PositionBoundaryCondition(NodeGrid,0);
    
    %     Reconfiguration of matrices
    % [ A_bb A_bf;
    %   A_fb A_ff]
    K= MatrixReconfiguration(K, boundNodes);
    M= MatrixReconfiguration(M, boundNodes);
    D= MatrixReconfiguration(D, boundNodes);
    
    % Extract submatrices
    M_FF = M((length(boundNodes)+1):end,(length(boundNodes)+1):end);
    M_FB = M((length(boundNodes)+1):end,1:length(boundNodes));

    D_FF = D((length(boundNodes)+1):end,(length(boundNodes)+1):end);
    D_FB = D((length(boundNodes)+1):end,1:length(boundNodes));
    
    K_FF = K((length(boundNodes)+1):end,(length(boundNodes)+1):end);
    K_FB = K((length(boundNodes)+1):end,1:length(boundNodes));
    invK_FF = inv(K_FF);
    % Eigenmodes of the free system

    [Phi, Omega] = eig(K_FF,M_FF,'vector');

    [Omega, ind] = sort(Omega);
    Phi = Phi(:, ind);
    Phi = [Phi(:,1:NumberOfModes), AdditionalModes];
    Omega = Phi.'*K_FF*Phi;
    % Calculate Matrices
    D_tilde = Phi.'*D_FF*Phi;

    M_bar = M_FF*(invK_FF\K_FB)-M_FB;
    D_bar = D_FF*(K_FF\K_FB)-D_FB;
    invK_FFK_FB = (K_FF\K_FB);



end