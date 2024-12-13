function [did_pass,too_many_inside,fig_handle] = NIS(innovation_vec,Sk,alpha,show_plot)
% Inputs:
% - total measurement state yhat(k) = ynom(k) + dyhat(k)
% - ytruth, yhat_plus = p x length of time array x N
% - N = number MC simulation runs
% - Sk = n x n x length of time array x N
% - alpha = scalar = significance level
% - show_plot = optional, if true then plot result
% 
% Outputs:
% - did_pass = epsilon_yk bar inside [r1,r2] for at least around
%       (1-alpha)*100% of time (this is good)
% - too_many_inside = = epsilon_yk bar inside [r1,r2] more than
%       (1-alpha)*100% of time (this is less good)
% - fig_handle = handle of figure if show_plot is true

p = size(innovation_vec,1); % n = number of states
ktot= size(innovation_vec,2); % length of time array
N = size(innovation_vec,3); % number of MC sims

% errors
e_yk = innovation_vec;

% calculate mean normalized estimation error squared
epsbar_yk = nan(1,ktot);
for ik = 1:ktot % for each time
    eps_yk = nan(1, N);
    for iMC = 1:N
        this_e_yk = squeeze( e_yk(:,ik,iMC) );
        thisSk = squeeze( Sk(:,:,ik,iMC) );
        eps_yk(iMC) = this_e_yk'*inv(thisSk)*this_e_yk;
    end
    epsbar_yk(ik) = 1/N*sum(eps_yk);
end

% find bounds for NEES test
r1 = chi2inv(alpha/2,N*p)/N;
r2 = chi2inv(1-alpha/2,N*p)/N;

inside_vals = sum(epsbar_yk>r1 & epsbar_yk<r2); % number of vals inside interval [r1,r2]
predicted_inside_vals = (1-alpha)*ktot; % number of vals we SHOULD see inside interval

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
    scatter(1:ktot, epsbar_yk,Color='b',DisplayName='NIS val')
    hold on
    plot([0, ktot],[r1 r1],'--',Color='r',DisplayName='r_1 bound')
    plot([0, ktot],[r2 r2],'--',Color='r',DisplayName='r_2 bound')
    hold off
    title('NIS Estimation Results')
    ylabel('NIS stat, $\bar{\epsilon_k}$',Interpreter='latex')
    xlabel('time step, k')
    legend
end

