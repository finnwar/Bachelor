function [mean_abs_error, mean_abs_squared_error, root_mean_squared_error, mean_rel_error,t] = StressErrorCalculation(t_ref,Sigma_ref,t_comp,Sigma_comp)
if length(t_ref) > length(t_comp)
Sigma_interp = zeros(length(Sigma_ref(:,1)),length(t_ref));
    for i= 1:length(Sigma_ref(:,1))
        Sigma_interp(i,:) = interp1(t_comp,Sigma_comp(i,:),t_ref,'spline');
        
    end
    t = t_ref;
    Sigma_comp = Sigma_interp;
else
Sigma_interp = zeros(length(Sigma_ref(:,1)),length(t_comp));
    for i = 1:length(Sigma_ref(:,1))
        Sigma_interp(i,:) = interp1(t_ref,Sigma_ref(i,:),t_comp,'spline');
    end
    t = t_comp;
    Sigma_ref = Sigma_interp;
end

abs_error = vecnorm(Sigma_comp-Sigma_ref);

rel_error = vecnorm(Sigma_comp-Sigma_ref)./vecnorm(Sigma_ref);

mean_abs_error = mean(abs_error);

mean_abs_squared_error = mean(abs_error.^2);

root_mean_squared_error = vecnorm(abs_error);

mean_rel_error = mean(rel_error,'all','omitmissing');

end