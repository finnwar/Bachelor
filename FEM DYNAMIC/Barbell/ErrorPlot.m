function ErrorPlot(U_ref,t_ref,U_comp,t_comp,Title)

if length(t_ref) > length(t_comp)
    for i= 1:length(U_ref(:,1))
        U_interp(i,:) = interp1(t_comp,U_comp(i,:),t_ref);
        
    end
    t = t_ref;
    U_comp = U_interp;
else
    for i = 1:length(U_ref(:,1))
        U_interp(i,:) = interp1(t_ref,U_ref(i,:),t_comp);
    end
    t = t_comp;
    U_ref = U_interp;
end

abs_error = abs(vecnorm(U_comp)-vecnorm(U_ref));

rel_error = vecnorm(U_comp)./vecnorm(U_ref);

% mean_abs_error = mean(abs_error);
% 
% mean_abs_squared_error = mean(abs_error.^2);
% 
% root_mean_squared_error = vecnorm(abs_error);
% 
% mean_rel_error = mean(rel_error);

figure
title(Title)
hold on
yyaxis left;
plot(t,abs_error)
yyaxis right;
plot(t,rel_error)
legend('Absolute Error','Relative Error')
hold off

end