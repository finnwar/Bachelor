function [abs_error,rel_error,total_error,mean_rel_error]= ErrorCalculation(t_ref,U_ref,t_comp,U_comp)
if length(t_ref) > length(t_comp)
    for i= 1:length(U_ref(:,1))
        U_interp(i,:) = interp1(t_comp,U_comp(i,:),t_ref,'spline');
        
    end
    t = t_ref;
    U_comp = U_interp;
elseif length(t_ref)==length(t_comp)
    t = t_ref;
else
    for i = 1:length(U_ref(:,1))
        U_interp(i,:) = interp1(t_ref,U_ref(i,:),t_comp,'spline');
    end
    t = t_comp;
    U_ref = U_interp;
end

abs_error = vecnorm(U_comp-U_ref);

rel_error = vecnorm(U_comp-U_ref)./vecnorm(U_ref);

total_error = vecnorm(abs_error);

mean_rel_error = mean(rel_error,'omitmissing');

end