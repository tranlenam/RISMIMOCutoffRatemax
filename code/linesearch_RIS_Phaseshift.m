function [mytheta_new,step_size_theta_new] = linesearch_RIS_Phaseshift(Hdirt,P,H1t,H2t,mytheta_current,x,...
    N0,step_size_theta_current,delta,rho)
% mytheta_new: new phase shift
% step_size_theta_new : new step size after line search procedure
step_size_theta_new = step_size_theta_current;

% current channel
H = Hdirt + H2t*diag(mytheta_current)*H1t; 

objprev = computeobjective(H*P*x,N0); % the current value

for iLineSearch =0:30 % maximum 30 steps
    % gradient step
    mytheta_new = mytheta_current-step_size_theta_new*grad_phase(Hdirt,P,H1t,H2t,mytheta_current,x,N0);
    % project onto feasible set
    mytheta_new = mytheta_new./abs(mytheta_new);
    % update channel with new phase shift
    Hnew = Hdirt + H2t*diag(mytheta_new)*H1t;
    % compute the value of the new phase shift
    objnew = computeobjective(Hnew*P*x,N0);
    if ((objprev - objnew) >= delta*(norm(mytheta_new-mytheta_current)^2)) % if better solution found, then break
        break;
    else % if not, then reduce step size by rho
        step_size_theta_new = step_size_theta_new*rho;
    end
end
end