function [K_tilde, M_tilde, D_tilde, f_tilde,Phi, eigenFrequencies] = ModalReduction(K,M,D,f,NumberOfModes, AdditionalModes)

    if nargin == 5
        %Calculate eigenmodes and frequencies
        [eigenVectors,eigenFrequencies] = eig(K,M,'vector');
        % Sorting eigenmodes via their frequencies from low to high
        [eigenFrequencies, index] = sort(eigenFrequencies);
        eigenVectors = eigenVectors(:,index);
    
    
    else
        % Calculate eigenmodes and frequencies
        [eigenVectors,eigenFrequencies] = eig(K,M,'vector');
        % Sorting eigenmodes via their frequencies from low to high
        [eigenFrequencies, index] = sort(eigenFrequencies);
        eigenVectors = eigenVectors(:,index);


        % Modal redcution with eigenmodes
        Phi = [eigenVectors(:,1:NumberOfModes), AdditionalModes];
    
        M_tilde = Phi.'*M*Phi;
        K_tilde = Phi.'*M*Phi;
        D_tilde = Phi.'*D*Phi;
        f_tilde = Phi.'*f;
    
    end




end