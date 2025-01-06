function [EigenVectorDerivative] =  NelsonEigenmodeSensitivity(A, dAdp, lambda, X, Y)

F = zeros(size(X));
V = zeros(size(X(:,1)));
CalcMatrix = zeros(size(X));
EigenVectorDerivative = zeros(size(X));


for i = 1:length(X(:,1))
    F(:,i) = X(:,i)*(Y(:,i).'*dAdp*X(:,i))-dAdp*X(:,i);
    
    [~,k] = max(abs(X(:,i)).*abs(Y(:,i)));


    
    CalcMatrix(1:(k-1),1:(k-1))=A(1:(k-1),1:(k-1))-lambda(i)*eye(k-1);
    CalcMatrix((k+1):end,(k+1):end)=A((k+1):end, (k+1):end)-lambda(i)*eye(length(V)-k);
    CalcMatrix(1:(k-1),(k+1):end) = A(1:(k-1),(k+1):end);
    CalcMatrix((k+1):end,1:(k-1)) = A((k+1):end,1:(k-1));
    CalcMatrix(k,k)=1;

    F(k) = 0;
    
    EigenVectorDerivative(:,i) = CalcMatrix\F;


end    


end