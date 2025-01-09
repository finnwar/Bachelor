function InitialPlottingFunction(PositionGrid)
    
    [I, J] = size(PositionGrid);

    PosMatrix = zeros(I*J, 2);

    for i=1:(I/2)
        for j = 1:J
            PosMatrix(i+j*I, 1) = PositionGrid(2*i-1, j);
            PosMatrix(i+j*I, 2) = PositionGrid(2*i, j);
        end
    end

    figure
    scatter(PosMatrix(:,1), PosMatrix(:,2),'filled')
    xlabel('x')
    ylabel('y')
end