% 5044 Final Project
% Sarah Luettgen, Abby Rindfuss, and Lisa Ventura
% Cooperative Location

% Housekeeping
clear; 
clc
close all

%% Parts I, Problems 1 and 2.

% ---- nominal conditions -----

L = .5;
%phi_g is between -5*pi/12 to 5*pi/12
%v_gmax = 3;
%omega_g is between -pi/6 to pi/6
%v_a is between 10 and 20
xi_g0 = 10;
eta_g0 = 0;
theta_g0 = pi/2;
v_g0 = 2;  % UGV nominal speed (m/s)
phi_g0 = -pi/18;  % UGV nominal steering angle (rad)

xi_a0 = -60;
eta_a0 = 0;
theta_a0 = -pi/2;
v_a0 = 12; % UAV nominal s peed (m/s)
omega_a0 = pi/25;  % UAV nominal turning rate (rad/s)

% delta T sampling rate
dt = .1;
tf = 100; %seconds
tarr = 0:dt:(tf); % t vector

% nominal u(t) vector components 
u1 = v_g0;
u2 = phi_g0;
u3 = v_a0;
u4 = omega_a0;
unom_t0 = [u1;u2;u3;u4];

% nominal x(t) vector components 
x1 = xi_g0;
x2 = eta_g0;
x3 = theta_g0;
x4 = xi_a0;
x5 = eta_a0;
x6 = theta_a0;
xnom_t0 = [x1;x2;x3;x4;x5;x6];

% nominal solution
thetag_dot_nom = v_g0/L*tan(phi_g0);
theta_g_nom = theta_g0+thetag_dot_nom*tarr;
thetaa_dot_nom = omega_a0;
theta_a_nom = theta_a0+thetaa_dot_nom*tarr;

xnom = [xi_g0+v_g0/thetag_dot_nom*sin(theta_g_nom)-v_g0/thetag_dot_nom*sin(theta_g0);...
        eta_g0-v_g0/thetag_dot_nom*cos(theta_g_nom)+v_g0/thetag_dot_nom*cos(theta_g0);...
        theta_g_nom;...
        xi_a0+v_a0/thetaa_dot_nom*sin(theta_a_nom)-v_a0/thetaa_dot_nom*sin(theta_a0);...
        eta_a0-v_a0/thetaa_dot_nom*cos(theta_a_nom)+v_a0/thetaa_dot_nom*cos(theta_a0);...
        theta_a_nom];

unom = repmat(unom_t0,1,length(tarr));

[Abar,Bbar,Cbar,Dbar,F,G,H,M] = get_dynamics_matrices(xnom,unom,dt,L);


%% Part I, Problem 3.

du = zeros(4,length(tarr)); % du vector
deltx0 = [0; 1; 0; 0; 0; 0.1]; %dx0
% deltx0 = [x1; x2; x3; x4; x5; x6];
d_state = nan(6,length(tarr)); %dx
d_state(:,1) = deltx0;

% solve for dx_k+1
for i = 2:length(tarr)
    F_t = squeeze(F(:,:,i));
    G_t = squeeze(G(:,:,i));
    d_state(:,i) = F_t*d_state(:,i-1) + G_t*du(:,i-1);
end

full_state = xnom + d_state; % x = x_nom + dx
vtilde_nom = zeros(size(full_state));
ynom = calc_obs_from_state(xnom,vtilde_nom);

% full state solved with ODE45
my_ode = @(t,y) NL_ode(t,y,v_g0,phi_g0,v_a0,omega_a0,[0;0;0],[0;0;0],L);
[t,xarr] = ode45(my_ode,tarr,xnom_t0+deltx0);
state_perturbed_ode=xarr';

[t,xarr] = ode45(my_ode,tarr,xnom_t0);
state_nom_ode=xarr';

% measurements
d_y = nan(5,length(tarr));
for i=1:length(tarr)
    M_t = squeeze(M(:,:,i));
    H_t = squeeze(H(:,:,i));
    d_y(:,i) = H_t*d_state(:,i)+M_t*du(:,i);
end

y_linearized = ynom+d_y;
y_ode = calc_obs_from_state(state_perturbed_ode,vtilde_nom);

% plotting
xunits = {'$\xi_g$ (m)','$\eta_g$ (m)','$\theta_g$ (rad)','$\xi_a$ (m)','$\eta_a$ (m)','$\theta_a$ (rad)'};
dxunits = {'$\delta \xi_g$ (m)','$\delta \eta_g$ (m)','$\delta \theta_g$ (rad)','$\delta \xi_a$ (m)','$\delta \eta_a$ (m)','$\delta \theta_a$ (rad)'};
yunits = {'$\gamma_{ag}$ (rad)','$\rho_{ga}$ (m)','$\gamma_{ga}$ (rad)','$\xi_a$ (m)','$\eta_a$ (m)'};
dyunits = {'$\delta \gamma_{ag}$ (rad)','$\delta \rho_{ga}$ (m)','$\delta \gamma_{ga}$ (rad)','$\delta \xi_a$ (m)','$\delta \eta_a$ (m)'};
wrap_indices_x = [3,6];
wrap_indices_y = [1,3];
% plot linearized perturbations
%fig1 = figure('units','normalized','outerposition',[0 1 .5 1]);
figure()
plot_states(tarr,d_state,dxunits,[])
sgtitle('Linearized Approximate Perturbations vs. Time','FontSize',14, 'Interpreter','latex')

% plot full linearized state
% fig2 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,full_state,xunits,wrap_indices_x)
sgtitle('States vs. Time, Linearized Approximate Dynamics Simulation','FontSize',14, 'Interpreter','latex')

% States vs. Time, Full Nonlinear Dynamics Simulation
% fig3 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,state_perturbed_ode,xunits,wrap_indices_x)
sgtitle('States vs. Time, Full Nonlinear Dynamics Simulation','FontSize',14, 'Interpreter','latex')

% plot measurements, linearized
%fig4 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,y_linearized,yunits,[1,3])
sgtitle('Linearized Approximate Dynamics Measurements','FontSize',14, 'Interpreter','latex')

% plot measurements from ode
%fig5 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,y_ode,yunits,[1,3])
sgtitle('Full Nonlinear Measurements','FontSize',14, 'Interpreter','latex')

% Full Nonlinear Dynamics Simulation minus full linearized state
%fig6 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,state_perturbed_ode-full_state,xunits,wrap_indices_x)
sgtitle('States vs. Time, ODE Minus Linearization','FontSize',14, 'Interpreter','latex')

% plot measurements from ode minus measurements, linearized
%fig7 = figure('units','normalized','outerposition',[0 1 .5 1]);
%   ^ This line of code isn't compatible with matlab 2024a, had to change it!
figure()
plot_states(tarr,y_ode-y_linearized,dyunits,[1,3])
sgtitle('ODE Measurements Minus Linearization Measurements','FontSize',14, 'Interpreter','latex')


%% Part II, Problem 4. 

% load the data, R, and Q matrices 
coopData = load('cooplocalization_finalproj_KFdata.mat');
Q = coopData.Qtrue;
R = coopData.Rtrue;
ydata = coopData.ydata;

MC_num = 100; % number of monte carlo simulations

% initialize covariance matrix 
    % 6x6, for each timestep, for each MC
P0 = eye(6); % initial state covariance matrix (IDK WHAT TO PUT HERE SO I MADE IT IDENTITY)
Pk_plus_all = zeros(6,6,length(tarr),MC_num); % Pk plus
Sk_plus_all = zeros(5,5,length(tarr),MC_num); % Pk plus

% initialize state matrix
    % 6x1, for each timestep, for each MC
dxhat_all = zeros(6,length(tarr),MC_num); % dx hat plus, given by KF
dx_truth_sim = zeros(6,length(tarr),MC_num); % simulated noisy truth perturbation
x_truth_sim = zeros(6,length(tarr),MC_num); % simulated truth state
xhat_all = zeros(6,length(tarr),MC_num);

% initialize measurements matrix
    % 5x1, for each timestep, for each MC
dy_KF = zeros(5, length(tarr),MC_num); % dy given by KF
dy_truth_sim = zeros(5, length(tarr),MC_num); % simulated measurement perturbations
y_truth_sim = zeros(5, length(tarr),MC_num); % ynom + noisy_dy
y_all = zeros(5, length(tarr),MC_num); % y given by KF

for m = 1:MC_num % for each MC iteration
    % ----------------
    % First, simulate truth state as "truth" for the NEES  test.
    % Also simulate corresponding measurements to use as "truth" for NIS 
    % test and as input to KF filter.
    xtrue0 = mvnrnd(xnom_t0,P0)'; % sample initial state % not sure if this is right!!
    dx_truth_sim(:,1,m) = xtrue0-xnom_t0;
    for k=2:length(tarr)
      
        %%simulate process noise and add to actual state
        wk = mvnrnd(zeros(1,6),Q)';
        F_t = squeeze(F(:,:,k));
        G_t = squeeze(G(:,:,k));

        dx_truth_sim(:,k,m) = F_t*dx_truth_sim(:,k-1,m) + G_t*du(:,k-1) + wk; 
        
        %%simulate measurement noise and add to sensor data
        vk = mvnrnd(zeros(1,size(ydata,1)),R)';
        H_t = squeeze(H(:,:,k));
        M_t = squeeze(M(:,:,k));
        dy_truth_sim(:,k,m) = H_t*dx_truth_sim(:,k,m) + M_t*du(:,k) + vk; 
    end
     % add simulated dx and dy to the nominal states
     x_truth_sim(:,:,m) = xnom + dx_truth_sim(:,:,m);
     y_truth_sim(:,:,m) = ynom + dy_truth_sim(:,:,m);
    
    % ----------------
    % Next, use the simulated "truth" measurements as inputs to KF

    % initialize KF
    dxhat0 = deltx0; % initial perturbation state estimate
    dxhat = zeros(6,length(tarr)); 
    dxhat(:,1) = dxhat0;

    dyhat = zeros(5,length(tarr));

    du0 = [0;0;0;0]; % fill in with real numbers later!!
    du = zeros(4,length(tarr));
    du(:,1) = du0;

    P0 = eye(6); % initial state covariance matrix (IDK WHAT TO PUT HERE SO I MADE IT IDENTITY)
    Pk = zeros(6,6,length(tarr));
    Pk(:,:,1) = P0;

    for k = 2:length(tarr) % for each timestep k 
        % must re-calculate F, G, H, Omega at each timestep! These are not
        % the same as the F, G, H Jacobians calculated in part 1. 
        Ftild_k = eye(6) + dt*Abar(:,:,k-1);
        Htild_k = H(:,:,k); % H tilde = H bar = C bar
        Gtild_k = dt*Bbar(:,:,k-1);
        Gamma = eye(6);
        Omegatild_k = dt*Gamma;

        % prediction step
        dxhat_minus = Ftild_k*dxhat(:,k-1) + Gtild_k*du(:,k-1);
        Pk_minus = Ftild_k*Pk(:,:,k-1)*Ftild_k' + Omegatild_k*Q*Omegatild_k';
        % du(:,k) =  
        % I'm confused about finding du
        % I think we don't need to?
        
        % gain K
        Skval = Htild_k*Pk_minus*Htild_k' + R;
        Skval = 0.5*(Skval + Skval');

        K = Pk_minus*Htild_k' * inv(Htild_k*Pk_minus*Htild_k' + R);

        % correction step
        Pk_plus = (eye(6) - K*Htild_k)*Pk_minus;

        pretend_y_data = y_truth_sim(:,k,m);
        predicted_y = ynom(:,k);
        %dy = ydata(:,k) - y_truth_sim(:,k,m);
        dy= pretend_y_data-predicted_y;
        innovation_k = dy-Htild_k*dxhat_minus;

        dxhat_plus = dxhat_minus + K*(innovation_k);

        % results of the state & covariance
        dxhat(:,k) = dxhat_plus;
        Pk(:,:,k) = Pk_plus;
        dyhat(:,k) = dy;
        Sk(:,:,k) = Skval;
        innovation(:,k) = innovation_k;

    end
    dxhat_all(:,:,m) = dxhat;
    Pk_plus_all(:,:,:,m) = Pk;
    dy_KF(:,:,m) = dyhat;
    xhat_all(:,:,m) = dxhat + x_truth_sim(:,:,m);
    y_all(:,:,m) = dyhat + y_truth_sim(:,:,m);
    Sk_all(:,:,:,m) = Sk;
    innovation_all(:,:,m) = innovation;
end

%% Plots for Problem 4a
% Plots for a single ‘typical’ simulation instance, showing the noisy simulated ground truth
% states, noisy simulated data, and resulting linearized KF state estimation errors
    % just picking monte carlo iteration #5 arbitrarily as the one to plot
% noisy simulated ground truth states + corresponding KF estimation 
figure()
plot_KF(tarr,x_truth_sim(:,:,5), dxhat_all(:,:,5), xunits, wrap_indices_x)
sgtitle('Simulated States, Linearized KF','FontSize',14, 'Interpreter','latex')

% noisy simulated data + corresponding KF estimation
figure()
plot_KF(tarr, y_truth_sim(:,:,5), y_all(:,:,5), yunits, wrap_indices_y)
sgtitle('Simulated Measurements, Linearized KF','FontSize',14, 'Interpreter','latex')


%% NEES test
xhat_plus = repmat(xnom,1,1,MC_num) + dxhat_all;
alpha_NEES = 0.05;
[did_pass_NEES,too_many_inside_NEES,fig_handle_NEES] = NEES(x_truth_sim, xhat_plus,Pk_plus_all,alpha_NEES,1);
% Inputs:
% - total state xhat^plus(k) = xnom(k) + dxhat^plus(k)
% - xtruth, xhat_plus = n x length of time array x N
% - N = number MC simulation runs
% - Pk_plus = n x n x length of time array x N
% - alpha = scalar = significance level
% - show_plot = optional, if true then plot result

%% NIS test
alpha_NIS = alpha_NEES;
[did_pass_NIS,too_many_inside_NIS,fig_handle_NIS] = NIS(innovation_all,Sk_all,alpha_NIS,1);
% Inputs:
% - total measurement state yhat(k) = ynom(k) + dyhat(k)
% - ytruth, yhat_plus = p x length of time array x N
% - N = number MC simulation runs
% - Sk = n x n x length of time array x N
% - alpha = scalar = significance level
% - show_plot = optional, if true then plot result


%% validation
for i=1:length(tarr)
    good_times = [.1,12.4,24.3];
    if not(isempty(find(good_times==tarr(i))))
        for j = 1:size(y_linearized,1)
            if (j==3 || j==1)
                toplot1 = mod(y_linearized(j,i)+pi,2*pi)-pi;
                toplot2 = mod(y_ode(j,i)+pi,2*pi)-pi;
            else
                toplot1= y_linearized(j,i);
                toplot2 = y_ode(j,i);
            end
            sprintf('%d, %2.2f, %2.2f, %2.2f',j,tarr(i),toplot1,toplot2)
        end
    end
end
%% Functions

% -------- ODE45 ----------

function yd = NL_ode(t,y,vg,phi,va,wa,w_tild_g,w_tild_a,L)
% full ODE for the dynamical system to be used in ODE 45
    xi_g=y(1);
    etag=y(2);
    theta_g=y(3);
    xi_a=y(4);
    etaa=y(5);
    theta_a=y(6);
    
    w_tild_xg = w_tild_g(1);
    w_tild_yg = w_tild_g(2);
    w_tild_wg = w_tild_g(3);

    w_tild_xa = w_tild_a(1);
    w_tild_ya = w_tild_a(2);
    w_tild_wa = w_tild_a(3);


    xi_g_dot = vg*cos(theta_g)+w_tild_xg;
    etag_dot = vg*sin(theta_g)+w_tild_yg;
    theta_g_dot = vg/L*tan(phi)+w_tild_wg;

    xi_a_dot = va*cos(theta_a)+w_tild_xa;
    etaa_dot = va*sin(theta_a)+w_tild_ya;
    theta_a_dot = wa+w_tild_wa;
    
    yd = [xi_g_dot; etag_dot; theta_g_dot; xi_a_dot; etaa_dot; theta_a_dot];
end

% --------- CT and DT MATRICES --------

function [Abar,Bbar,Cbar,Dbar,F,G,H,M] = get_dynamics_matrices(xnom,unom,dt,L)
% get matrices for the dynamical system
% xnom = xnom(t) where xnom is size 6x[length of time array]
% unom = unom(t) where unom is size 4x[length of time array]
% matrixes are size rows x cols x [length of time array]

num_times = size(xnom,2);
F = nan(6,6,num_times);
G = nan(6,4,num_times);
H = nan(5,6,num_times);
M = nan(5,4,num_times);
Abar = nan(6,6,num_times);
Bbar = nan(6,4,num_times);
Cbar = nan(5,6,num_times);
Dbar = nan(5,4,num_times);

    % ------- Jacobians -------
for ti = 1:num_times
    x1 = xnom(1,ti);
    x2 = xnom(2,ti);
    x3 = xnom(3,ti);
    x4 = xnom(4,ti);
    x5 = xnom(5,ti);
    x6 = xnom(6,ti);
    u1 = unom(1,ti);
    u2 = unom(2,ti);
    u3 = unom(3,ti);
    u4 = unom(4,ti);
    Abar(:,:,ti) = [0 0 -u1*sin(x3) 0 0 0; ...
            0 0 u1*cos(x3) 0 0 0;...
            0 0 0 0 0 0; ...
            0 0 0 0 0 -u3*sin(x6); ...
            0 0 0 0 0 u3*cos(x6); ...
            0 0 0 0 0 0];
    
    
    Bbar(:,:,ti) = [cos(x3) 0 0 0; ...
            sin(x3) 0 0 0; ...
            (1/L)*tan(u2) u1/L*(sec(u2))^2 0 0; ...
            0 0 cos(x6) 0; ...
            0 0 sin(x6) 0; ...
            0 0 0 1];
    
    abv = (x4-x1)^2 + (x5-x2)^2;
    
    Cbar(:,:,ti) = [(x5-x2)/abv (x1-x4)/abv -1 (x2-x5)/abv (x4-x1)/abv 0; ...
            (x1-x4)/sqrt(abv) (x2-x5)/sqrt(abv) 0 (x4-x1)/sqrt(abv) (x5-x2)/sqrt(abv) 0; ...
            (x5-x2)/abv (x1-x4)/abv 0 (x2-x5)/abv (x4-x1)/abv 0; ...
            0 0 0 1 0 0; ...
            0 0 0 0 1 0];


    Dbar(:,:,ti) = zeros(5,4);
    
    
    % Part I, Problem 2.
    
    z = [Abar(:,:,ti) Bbar(:,:,ti); zeros(4,6) zeros(4)];
    ez = expm(z*dt);
    
    F(:,:,ti) = ez(1:6, 1:6);
    G(:,:,ti) = ez(1:6, 7:10);
    H(:,:,ti) = Cbar(:,:,ti);
    M(:,:,ti) = Dbar(:,:,ti);
end

end

% ------ PLOTTING FUNCTION -------

function plot_states(tarr,state,ylabels,wrap_indices)
% plot the states using subplots
    for iw = 1:length(wrap_indices)
        state(wrap_indices(iw),:)= mod(state(wrap_indices(iw),:)+pi,2*pi)-pi;  
    end

    for i=1:size(state,1)
        subplot(size(state,1),1,i)
        plot(tarr,state(i,:),'Color','blue','LineWidth',1.5)
        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
        grid on
    end

end


% ------ KALMAN FILTER PLOTTING FUNCTION -------

function plot_KF(tarr,sim_state,KF_state,ylabels,wrap_indices)
% plot the states using subplots
    for iw = 1:length(wrap_indices)
        sim_state(wrap_indices(iw),:)= mod(sim_state(wrap_indices(iw),:)+pi,2*pi)-pi;  
        KF_state(wrap_indices(iw),:)= mod(KF_state(wrap_indices(iw),:)+pi,2*pi)-pi; 
    end

    for i=1:size(sim_state,1)
        subplot(size(sim_state,1),1,i)
        plot(tarr,sim_state(i,:),'Color','blue','LineWidth',1.5)
        hold on
        plot(tarr,KF_state(i,:),'Color','red','LineWidth',1.5)
        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
        grid on
        legend('Simulated ground truth', 'KF output')
    end

end

function y = calc_obs_from_state(x,vtilde)
    x1 = x(1,:)+vtilde(1,:);
    x2 = x(2,:)+vtilde(1,:);
    x3 = x(3,:)+vtilde(1,:);
    x4 = x(4,:)+vtilde(1,:);
    x5 = x(5,:)+vtilde(1,:);
    x6 = x(6,:)+vtilde(1,:);
    %x3= mod(x3+pi,2*pi)-pi;
    %x6= mod(x6+pi,2*pi)-pi;
       
    y = [atan2( (x5-x2) , (x4-x1) ) - x3;...
             sqrt( (-x4+x1).^2 + (-x5+x2).^2 );...
             atan2( (x2-x5) , (x1-x4) ) - x6;...
             x4;...
             x5];

end