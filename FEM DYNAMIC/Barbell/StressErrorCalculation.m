function [mean_abs_error, mean_abs_squared_error, root_mean_squared_error, mean_rel_error] = ErrorCalculation(t_ref,Sigma_ref,t_comp,Sigma_comp)
if length(t_ref) > length(t_comp)
    for i= 1:length(Sigma_ref(:,1))
        Sigma_interp(i,:) = interp1(t_comp,Sigma_comp(i,:),t_ref);
        
    end
    t = t_ref;
    Sigma_comp = Sigma_interp;
else
    for i = 1:length(Sigma_ref(:,1))
        Sigma_interp(i,:) = interp1(t_ref,Sigma_ref(i,:),t_comp);
    end
    t = t_comp;
    Sigma_ref = Sigma_interp;
end

abs_error = vecnorm(Sigma_comp-Sigma_ref);

rel_error = abs((Sigma_comp-Sigma_ref)./Sigma_ref);

mean_abs_error = mean(abs_error);

mean_abs_squared_error = mean(abs_error.^2);

root_mean_squared_error = vecnorm(abs_error);

mean_rel_error = mean(rel_error,'all','omitmissing');

end