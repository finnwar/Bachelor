function PlotDisplacement(U,NodeGrid,NodePosition)

    Position = zeros(length(U), 1);

    for i = 1:length(NodeGrid(1,:))
        for j = 1:length(NodeGrid(:,1))
            Position(NodeGrid(j,i))=U(NodeGrid(j,i))+NodePosition(j,i);
        end
    end
PositionXY = zeros(length(U)/2, 2);
PositionXY(1:end, 1) = Position(1:2:(end-1));
PositionXY(1:end, 2) = Position(2:2:end);

    figure
    scatter(PositionXY(:,1),PositionXY(:,2),'filled')

end
