function B = B_matrix(xi, eta, dxidx, dxidy, detadx, detady)
    
    [~,dNdxi,dNdeta] = ShapeFunctions(xi,eta);
    Temp = [dxidx detadx; dxidy detady]*[dNdxi;dNdeta];
    
    dNdx=Temp(1,:);
    dNdy=Temp(2,:);

    B = [dNdx(1) 0 dNdx(2) 0 dNdx(3) 0 dNdx(4) 0; ...
        0 dNdy(1) 0 dNdy(2) 0 dNdy(3) 0 dNdy(4); ...
         dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];
end