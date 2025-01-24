function A_tilde = MatrixReconfiguration(A, boundaryNodes)
A_tilde = zeros(size(A));

A_tilde(1:length(boundaryNodes),:) = A(boundaryNodes,:);
A_tilde(:,1:length(boundaryNodes)) = A(:,boundaryNodes);
A_tilde((length(boundaryNodes)+1):end,:) = A(setdiff(1:length(A(1,:)),boundaryNodes),:);
A_tilde(:,(length(boundaryNodes)+1):end) = A(:,setdiff(1:length(A(1,:)),boundaryNodes));
end