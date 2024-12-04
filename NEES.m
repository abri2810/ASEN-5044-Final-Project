function [did_pass,too_many_inside,fig_handle] = NEES(xtruth, xhat_plus,Pk_plus,alpha,show_plot)
% Inputs:
% - total state xhat^plus(k) = xnom(k) + dxhat^plus(k)
% - xtruth, xhat_plus = n x length of time array x N
% - N = number MC simulation runs
% - Pk_plus = n x n x length of time array x N
% - alpha = scalar = significance level
% - show_plot = optional, if true then plot result
% 
% Outputs:
% - did_pass = epsilon_xk bar inside [r1,r2] for around
%       (1-alpha)*100% of time (this is good)
% - too_many_inside = = epsilon_xk bar inside [r1,r2] more than
%       (1-alpha)*100% of time (this is less good)
% - fig_handle = handle of figure if show_plot is true

nstates = size(xhat_plus,1); % n = number of states
ktot= size(xhat_plus,2); % length of time array
N = size(xhat_plus,3); % number of MC sims

% errors
e_xk = xtruth-xhat_plus;

% calculate mean normalized estimation error squared
epsbar_xk = nan(1,ktot);
for ik = 1:ktot % for each time
    eps_xk = nan(1,N);
    for iMC = 1:N
        this_e_xk = squeeze( e_xk(:,ik,iMC) );
        thisPk_plus = squeeze( Pk_plus(:,:,ik,iMC) );
        eps_xk(iMC) = this_e_xk'*inv(thisPk_plus)*this_e_xk;             
    end
    epsbar_xk(ik) = mean(eps_xk);
end

% find bounds for NEES test
r1 = chi2inv(alpha/2,N*nstates)/N;
r2 = chi2inv(1-alpha/2,N*nstates)/N;

inside_vals = sum(epsbar_xk>r1 & epsbar_xk<r2); % number of vals inside interval [r1,r2]
predicted_inside_vals = (1-alpha)*100*ktot; % number of vals we SHOULD see inside interval

% test to see if KF filter passes NEES test
% may need to make this first "if" statement more lenient 
if (inside_vals<ceil(predicted_inside_vals) && inside_vals>floor(predicted_inside_vals) )
    % about (1-alpha)*100% of epsilon bar values are inside expected bounds
    did_pass = 1; 
    too_many_inside = 0;
elseif inside_vals>ceil(predicted_inside_vals)
    % > (1-alpha)*100% of epsilon bar values are inside expected bounds
    did_pass = 0; 
    too_many_inside = 1;
else
    % < (1-alpha)*100% of epsilon bar values are inside expected bounds
    did_pass = 0;
    too_many_inside = 0;
end

if not(exist('show_plot','var')) show_plot = 0; end

if show_plot
    fig_handle = figure();
    scatter(1:ktot, epsbar_xk,Color='b',DisplayName='NEES val')
    hold on
    plot([0, ktot],[r1 r1],'--',Color='r',DisplayName='r_1 bound')
    plot([0, ktot],[r2 r2],'--',Color='r',DisplayName='r_2 bound')
    hold off
    title('NEES Estimation Results')
    ylabel('NEES stat, $\bar{\epsilon_k}$',Interpreter='latex')
    xlabel('time step, k')
    legend
end

