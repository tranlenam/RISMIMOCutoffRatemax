function [P_new,step_size_P_new] = linesearch_Tx_Precoder(Hdirt,P_current,H1t,H2t,RIS_phase,x,...
    N0,step_size_P_current,delta,rho,Pow)
% RIS_phase_new: new phase shift
% step_size_theta_new : new step size after line search procedure
step_size_P_new = step_size_P_current;

% current channel
H = Hdirt + H2t*diag(RIS_phase)*H1t; 
objprev = computeobjective(H*P_current*x,N0); % the current objective value

for iLineSearch =0:30 % maximum 30 steps
    % gradient step
    P_new = P_current-step_size_P_new*grad_P(Hdirt,P_current,H1t,H2t,RIS_phase,x,N0);
    % project onto feasible set
    if (norm(P_new,'fro')>sqrt(Pow))
        P_new = P_new/norm(P_new)*sqrt(Pow);
    end
    % compute the value of the new phase shift
    objnew = computeobjective(H*P_new*x,N0);
    if ((objprev - objnew) >= delta*(norm(P_new-P_current,'fro')^2)) % if better solution found, then break
        break;
    else % if not, then reduce step size by rho
        step_size_P_new = step_size_P_new*rho;
    end
end
end