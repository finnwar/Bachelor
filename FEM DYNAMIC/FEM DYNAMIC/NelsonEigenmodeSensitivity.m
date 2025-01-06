function [EigenVectorDerivative] =  NelsonEigenmodeSensitivity(A, dAdp, lambda, X, Y)

F = zeros(size(X));
V = zeros(size(X(:,1)));
CalcMatrix = zeros(size(X));
EigenVectorDerivative = zeros(size(X));

M = eye(size(X));           % Still have to figure out how to get M and dMdp
dMdp = eye(size(X));

for i = 1:length(X(:,1))
    F(:,i) = X(:,i)*(Y(:,i).'*dAdp*X(:,i))-dAdp*X(:,i);
    
    [~,k] = max(abs(X(:,i)).*abs(Y(:,i)));          %if chosen x_k << max(x) and y_k<<max(y)
    

    
    CalcMatrix(1:(k-1),1:(k-1))=A(1:(k-1),1:(k-1))-lambda(i)*eye(k-1);
    CalcMatrix((k+1):end,(k+1):end)=A((k+1):end, (k+1):end)-lambda(i)*eye(length(V)-k);
    CalcMatrix(1:(k-1),(k+1):end) = A(1:(k-1),(k+1):end);
    CalcMatrix((k+1):end,1:(k-1)) = A((k+1):end,1:(k-1));
    CalcMatrix(k,k)=1;

    F(k) = 0;

    V = CalcMatrix\F;
    c = -real(X(:,i)'*M*V)+0.5*X(:,i)'*dMdp*X(:,i);
    EigenVectorDerivative(:,i) = V + c*X(:,i); 

end    


end