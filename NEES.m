function [did_pass,too_many_inside] = NEES(xtruth, xhat_plus,Pk_plus,alpha)
% Inputs:
% - total state xhat^plus(k) = xnom(k) + dxhat^plus(k)
% - xtruth, xhat_plus = n x length of time array x N
% - N = number MC simulation runs
% - Pk_plus = n x n x length of time array x N
% - alpha = scalar = significance level
% 
% Outputs:
% - did_pass = epsilon_xk bar inside [r1,r2] for at least around
%       (1-alpha)*100% of time (this is good)
% - too_many_inside = = epsilon_xk bar inside [r1,r2] more than
%       (1-alpha)*100% of time (this is less good)

nstates = size(xhat_plus,1); % n = number of states
ktot= size(xhat_plus,2); % length of time array
N = size(xhat_plus,3); % number of MC sims

% errors
e_xk = xtruth-xhat_plus;

% calculate mean normalized estimation error squared
epsbar_xk = nan(1,ktot);
for ik = 1:ktot % for each time
    eps_xk = nan(nstates, N);
    for iMC = 1:N
        this_e_xk = squeeze( e_xk(:,ik,iMC) );
        thisPk_plus = squeeze( Pk_plus(:,:,ik,iMC) );
        eps_xk(iMC) = this_e_xk'*inv(thisPk_plus)*this_e_xk;
    end
    epsbar_xk(ik) = 1/N*sum(eps_xk);
end

r1 = chi2inv(alpha/2,N*nstates)/N;
r2 = chi2inv(1-alpha/2,N*nstates)/N;
inside_vals = sum(epsbar_xk>r1 & epsbar_xk<r2); % number of vals inside interval [r1,r2]
predicted_inside_vals = (1-alpha)*100*ktot; % number of vals we SHOULD see inside interval
if (inside_vals<ceil(predicted_inside_vals) && inside_vals>floor(predicted_inside_vals) )
    did_pass = 1;
    too_many_inside = 0;
elseif inside_vals>ceil(predicted_inside_vals)
    did_pass = 1; 
    too_many_inside = 1;
else
    did_pass = 0;
    too_many_inside = 0;
end


