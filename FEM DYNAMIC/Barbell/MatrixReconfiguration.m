function [A_tilde, A_ee, A_ei, A_ie, A_ii] = MatrixReconfiguration(A, boundaryNodes)


A_ii= A;

A_ii(:,boundaryNodes) = [];
A_ii(boundaryNodes,:) = [];

A_ee = A;
A_ee(:,setdiff(1:length(A(1,:)),boundaryNodes)) = [];
A_ee(setdiff(1:length(A(1,:)),boundaryNodes),:) = [];

A_ei = A(boundaryNodes,setdiff(1:length(A(1,:)),boundaryNodes));

A_ie = A(setdiff(1:length(A(1,:)),boundaryNodes),boundaryNodes);

A_tilde = [A_ee A_ei;
           A_ie A_ii];

end