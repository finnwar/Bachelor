%Shape functions of linear 4-Node Qudrilateral

function [N, dNdxi, dNdeta]= ShapeFunctions(xi,eta)
   N = [(1/4)*(1-xi)*(1-eta), (1/4)*(1+xi)*(1-eta), (1/4)*(1+xi)*(1+eta), (1/4)*(1-xi)*(1+eta)];
   dNdxi = [(-1/4)*(1-eta), (1/4)*(1-eta), (1/4)*(1+eta),(-1/4)*(1+eta)];
   dNdeta = [(-1/4)*(1-xi), (-1/4)*(1+xi), (1/4)*(1+xi), (1/4)*(1-xi)];
end


