% 5044 Final Project
% Sarah Luettgen, Abby Rindfuss, and Lisa Ventura
% Cooperative Location

% Housekeeping
clear; 
clc
close all

rng(50)

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

% xnom = [xi_g0+v_g0/thetag_dot_nom*sin(theta_g_nom)-v_g0/thetag_dot_nom*sin(theta_g0);...
%         eta_g0-v_g0/thetag_dot_nom*cos(theta_g_nom)+v_g0/thetag_dot_nom*cos(theta_g0);...
%         theta_g_nom;...
%         xi_a0+v_a0/thetaa_dot_nom*sin(theta_a_nom)-v_a0/thetaa_dot_nom*sin(theta_a0);...
%         eta_a0-v_a0/thetaa_dot_nom*cos(theta_a_nom)+v_a0/thetaa_dot_nom*cos(theta_a0);...
%         theta_a_nom];

% Calculate xnom (an alternative to the above matrix)
my_ode = @(t,y) NL_ode(t,y,v_g0,phi_g0,v_a0,omega_a0,[0;0;0],[0;0;0],L);
[t,xarr] = ode45(my_ode,tarr,xnom_t0);
xnom=xarr';

unom = repmat(unom_t0,1,length(tarr));

[Abar,Bbar,Cbar,Dbar,F,G,H,M] = get_dynamics_matrices(xnom,unom,dt,L);


%% Part I, Problem 3.

du = zeros(4,length(tarr)); % du vector
deltx0 = [0; 1; 0; 0; 0; 0.1]; %dx(t=0)
d_state = nan(6,length(tarr)); %dx(t)
d_state(:,1) = deltx0;

% solve for dx_k+1
for i = 2:length(tarr)
    F_t = squeeze(F(:,:,i));
    G_t = squeeze(G(:,:,i));
    d_state(:,i) = F_t*d_state(:,i-1) + G_t*du(:,i-1);
end

full_state = xnom + d_state; % x = x_nom + dx
vtilde_nom = zeros(5,size(full_state,2)); % measurement noise
ynom = calc_obs_from_state(xnom,vtilde_nom); % ynom

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
figure(Visible="off")
plot_states(tarr,d_state,dxunits,[])
sgtitle('Linearized Approximate Perturbations vs. Time','FontSize',14, 'Interpreter','latex')

% plot full linearized state
figure(Visible="off")
plot_states(tarr,full_state,xunits,wrap_indices_x)
sgtitle('States vs. Time, Linearized Approximate Dynamics Simulation','FontSize',14, 'Interpreter','latex')

% States vs. Time, Full Nonlinear Dynamics Simulation
figure(Visible="off")
plot_states(tarr,state_perturbed_ode,xunits,wrap_indices_x)
sgtitle('States vs. Time, Full Nonlinear Dynamics Simulation','FontSize',14, 'Interpreter','latex')

% plot measurements, linearized
figure(Visible="off")
plot_states(tarr,y_linearized,yunits,[1,3])
sgtitle('Linearized Approximate Dynamics Measurements','FontSize',14, 'Interpreter','latex')

% plot measurements from ode
figure(Visible="off")
plot_states(tarr,y_ode,yunits,[1,3])
sgtitle('Full Nonlinear Measurements','FontSize',14, 'Interpreter','latex')

% Full Nonlinear Dynamics Simulation minus full linearized state
figure(Visible="off")
plot_states(tarr,state_perturbed_ode-full_state,xunits,wrap_indices_x)
sgtitle('Nonlinear States Minus Linearized States','FontSize',14, 'Interpreter','latex')

% plot measurements from ode minus measurements, linearized
figure(Visible="off")
plot_states(tarr,y_ode-y_linearized,dyunits,[1,3])
sgtitle('Nonlinear Measurements Minus Linearization Measurements','FontSize',14, 'Interpreter','latex')


% More plots!
% comparing linear states with nonlinear states (just plotting both
% together) 
figure(Visible="off")
plot_both_states(tarr,state_perturbed_ode,full_state,xunits,wrap_indices_x)
legend('Nonlinear','Linearized','FontSize',12,'Interpreter','latex')
sgtitle('Nonlinear and Linearized States','FontSize',14, 'Interpreter','latex')

figure(Visible="off")
plot_both_states(tarr,y_ode,y_linearized,yunits,[1,3])
legend('Nonlinear','Linearized','FontSize',12,'Interpreter','latex')
sgtitle('Nonlinear and Linearized Measurements','FontSize',14, 'Interpreter','latex')


%% Part II, Problem 4. 

% load the data, R, and Q matrices 
coopData = load('cooplocalization_finalproj_KFdata.mat');
Q = coopData.Qtrue;
R = coopData.Rtrue;

MC_num = 100; % number of monte carlo simulations

% initialize covariance matrix 
    % 6x6, for each timestep, for each MC
Pk_plus_all = zeros(6,6,length(tarr),MC_num); % Pk plus
Sk_all = zeros(5,5,length(tarr),MC_num); % Sk
innovation_all = zeros(5,length(tarr),MC_num); % innovation vextor

% initialize state matrix
    % 6x1, for each timestep, for each MC
dxhat_all = zeros(6,length(tarr),MC_num); % dx hat plus, given by KF
%dx_truth_sim = zeros(6,length(tarr),MC_num); % simulated noisy truth perturbation
x_truth_sim = zeros(6,length(tarr),MC_num); % simulated truth state
xhat_all = zeros(6,length(tarr),MC_num);

% initialize measurements matrix
    % 5x1, for each timestep, for each MC
dy_KF = zeros(5, length(tarr),MC_num); % dy given by KF
% dy_truth_sim = zeros(5, length(tarr),MC_num); % simulated measurement perturbations
y_truth_sim = zeros(5, length(tarr),MC_num); % ynom + noisy_dy
y_all = zeros(5, length(tarr),MC_num); % y given by KF

% initialize error bounds
xsigmas_all = zeros(6,length(tarr), MC_num);
ysigmas_all = zeros(5,length(tarr), MC_num);

% TUNING THE Q MATRIX
    % manually adjusting based on error plots and NEES plots
Q_LKF = diag([1000, 1000, .3, 100, 100, 2]);
R_LKF = R;
P0_LKF = diag([.5 .5 .03 300 300 .3]);


for m = 1:MC_num % for each MC iteration
    % --------------- 
    % Simulate truth state for NEES and NIS tests
    xtrue0 = mvnrnd(xnom_t0, Q)'; % initial state
    x_truth_sim(:,1,m)= xtrue0;

    
    % ----------------
    % Next, use the simulated "truth" measurements as inputs to KF

    % initialize KF
    dxhat0_LKF = deltx0; % initial perturbation state estimate
    du0_LKF = [0;0;0;0]; % zeros
    for k = 2:length(tarr) % for each timestep k 

        t1 = tarr(k-1);
        t2 = tarr(k);

        wk = mvnrnd(zeros(1, 6), Q)'; 
        vk = mvnrnd(zeros(1, 5), R)'; % Measurement noise

        my_ode = @(t,y) NL_ode(t,y,v_g0,phi_g0,v_a0,omega_a0,wk(1:3),wk(4:6),L);
        [t,x] = ode45(my_ode,[t1 t2],x_truth_sim(:,k-1,m)); 
        x_truth_sim(:,k,m) = x(end,:)';

        y_truth_sim(:,k,m) = calc_obs_from_state(x_truth_sim(:,k,m),vk);

    end
    
    % perform LKF
    [Pk,dyhat,dxhat,innovation,Sk,xsigmas,ysigmas] = LKF(y_truth_sim(:,:,m),dxhat0_LKF,P0_LKF,du0_LKF,tarr,Q_LKF,R_LKF,Abar,Bbar,H,xnom);

    % Save important values
    xsigmas_all(:,:,m) = xsigmas;
    ysigmas_all(:,:,m) = ysigmas;
    dxhat_all(:,:,m) = dxhat;
    Pk_plus_all(:,:,:,m) = Pk;
    dy_KF(:,:,m) = dyhat;
    xhat_all(:,:,m) = dxhat + xnom;
    y_all(:,:,m) = dyhat + ynom;
    Sk_all(:,:,:,m) = Sk;
    innovation_all(:,:,m) = innovation;
end

%% Plots for Problem 4a
% Plots for a single ‘typical’ simulation instance, showing the noisy simulated ground truth
% states, noisy simulated data, and resulting linearized KF state estimation errors
    % just picking monte carlo iteration #5 arbitrarily as the one to plot
% noisy simulated ground truth states + corresponding KF estimation 
figure()
plot_KF(tarr, x_truth_sim(:,:,5), xhat_all(:,:,5), xsigmas_all(:,:,5), xunits, wrap_indices_x)
sgtitle('Simulated States, Linearized KF','FontSize',14, 'Interpreter','latex')

% noisy simulated data + corresponding KF estimation
figure()
plot_KF(tarr, y_truth_sim(:,:,5), y_all(:,:,5), ysigmas_all(:,:,5), yunits, wrap_indices_y)
sgtitle('Simulated Measurements, Linearized KF','FontSize',14, 'Interpreter','latex')

% Plotting state (x) errors
figure()
plot_errors(tarr,x_truth_sim(:,:,5) - xhat_all(:,:,5),xsigmas_all(:,:,5),xunits)
sgtitle('State Errors, Linearized KF','FontSize',14, 'Interpreter','latex')

% Plotting measurement (y) errors
figure()
plot_errors(tarr(2:end),innovation_all(:,2:end,5),ysigmas_all(:,2:end,5), yunits)
sgtitle('Measurement Errors, Linearized KF','FontSize',14, 'Interpreter','latex')

% NEES test
xhat_plus = repmat(xnom,1,1,MC_num) + dxhat_all;
alpha_NEES = 0.05;
[did_pass_NEES,too_many_inside_NEES,fig_handle_NEES] = NEES(x_truth_sim, xhat_plus,Pk_plus_all,alpha_NEES,1);

% NIS test
alpha_NIS = alpha_NEES;
[did_pass_NIS,too_many_inside_NIS,fig_handle_NIS] = NIS(innovation_all,Sk_all,alpha_NIS,1);

%% Part II, Problem 5

% Initialization
Pk_plus_all = zeros(6, 6, length(tarr), MC_num); % State covariance
Sk_all = zeros(5, 5, length(tarr), MC_num);
y_KF = zeros(5, length(tarr), MC_num); % Predicted measurements
xhat_all = zeros(6, length(tarr), MC_num); % Final state estimate
y_act = zeros(5, length(tarr), MC_num); % Actual measurements
sigmas_all = zeros(6, 6, length(tarr), MC_num);

x_truth_sim = zeros(6,length(tarr), MC_num);
y_truth_sim = zeros(5, length(tarr),MC_num);

xhat0_EKF = xnom(:,1);
Pk0_EKF = diag([300 300 3 300 300 3]); 
Q_EKF = Q;

% initialize error bounds
    xsigmas_all = zeros(6,length(tarr), MC_num);
    ysigmas_all = zeros(5,length(tarr), MC_num);

for m = 1:MC_num % Monte Carlo iterations

    % Simulate truth state for NEES and NIS tests
    xtrue0 = mvnrnd(xnom_t0, Q)'; % initial state 
    x_truth_sim(:,1,m)= xtrue0;

    % Simulate measurements
    for k = 2:length(tarr)

        t1 = tarr(k-1);
        t2 = tarr(k);

        wk = mvnrnd(zeros(1, 6), Q)'; 
        vk = mvnrnd(zeros(1, 5), R)'; % Measurement noise

        my_ode = @(t,y) NL_ode(t,y,v_g0,phi_g0,v_a0,omega_a0,wk(1:3),wk(4:6),L);
        [t,x] = ode45(my_ode,[t1 t2],x_truth_sim(:,k-1,m)); 
        x_truth_sim(:,k,m) = x(end,:)';

        y_truth_sim(:,k,m) = calc_obs_from_state(x_truth_sim(:,k,m),vk);

    end


    % Initialize EKF

    [Pk_all, yhat,xhat_p,innovation,Sk_collect,sigmas_collect,xsigmas,ysigmas] = EKF(y_truth_sim(:,:,m),xhat0_EKF,Pk0_EKF,unom,L,tarr,Q_EKF,R);
    
    % Save results for this Monte Carlo iteration
    Pk_plus_all(:, :, :, m) = Pk_all;
    y_all(:, :, m) = yhat;
    xhat_all(:, :, m) = xhat_p;
    innovation_all(:,:,m) = innovation;
    Sk_all(:,:,:,m) = Sk_collect;

    xsigmas_all(:,:,m) = xsigmas;
    ysigmas_all(:,:,m) = ysigmas;

    sigmas_all(:,:,:,m) = sigmas_collect;
    
end


%% Plots for Problem 5a/EKF
% Plots for a single ‘typical’ simulation instance, showing the noisy simulated ground truth
% states, noisy simulated data, and resulting linearized KF state estimation errors
    % just picking monte carlo iteration #5 arbitrarily as the one to plot
% noisy simulated ground truth states + corresponding KF estimation 
figure()
plot_EKF(tarr,x_truth_sim(:,:,5), xhat_all(:,:,5), xunits, wrap_indices_x)
sgtitle('Simulated States, EKF','FontSize',14, 'Interpreter','latex')

% noisy simulated data + corresponding KF estimation
figure()
plot_EKF(tarr(2:end), y_truth_sim(:,2:end,5), y_all(:,2:end,5), yunits, wrap_indices_y)
sgtitle('Simulated Measurements, EKF','FontSize',14, 'Interpreter','latex')

% plot errors
figure()
error1 = x_truth_sim(:,:,5)-xhat_all(:,:,5);
plot_error(tarr(2:end), error1(:,2:end),sigmas_all(:,:,2:end,:), xunits)

% Plot ground truth positions
figure();
hold on;
plot(x_truth_sim(1,:,1), x_truth_sim(2,:,1), 'b', x_truth_sim(4,:,1), x_truth_sim(5,:,1), 'r');
% Add initial & final points
plot(x_truth_sim(1,1,1), x_truth_sim(2,1,1), 'bo', x_truth_sim(4,1,1), x_truth_sim(5,1,1), 'ro');
plot(x_truth_sim(1,end,1), x_truth_sim(2,end,1), 'bx', x_truth_sim(4,end,1), x_truth_sim(5,end,1), 'rx');
hold off;
xlabel('$\xi$ (m)','FontSize',12, 'Interpreter', 'latex');
ylabel('$\eta$ (m)', 'FontSize',12,'Interpreter', 'latex');
legend('Ground Vehicle', 'Air Vehicle','FontSize',12, 'Interpreter', 'latex');
title('GT Ground and Air Vehicle Positions');

% NEES test for EKF
xhat_plus = xhat_all; 
alpha_NEES = 0.05;
[did_pass_NEES,too_many_inside_NEES,fig_handle_NEES] = NEES(x_truth_sim, xhat_plus,Pk_plus_all,alpha_NEES,1);

% NIS test for EKF
alpha_NIS = alpha_NEES;
[did_pass_NIS,too_many_inside_NIS,fig_handle_NIS] = NIS(innovation_all,Sk_all,alpha_NIS,1);


%% Part II number 6
% Initialization
ydata = coopData.ydata;
tvec = coopData.tvec;

% EKF
[Pk_all, yhat,xhat,innovation,Sk_collect,sigmas_collect,xsigmas,ysigmas] = EKF(ydata,xhat0_EKF,Pk0_EKF,unom,L,tvec,Q_EKF,R);

% measurements + corresponding KF estimation
fig_EKF_yadata = figure();
plot_EKF(tvec(2:end), ydata(:,2:end), yhat(:,2:end), yunits, wrap_indices_y,1)
sgtitle('Measurements, EKF for ydata','FontSize',14, 'Interpreter','latex')
saveas(fig_EKF_yadata, 'images/EKF_ydata_pt6.png')

% plot 2*sigma_x
xsigmas_EKF = xsigmas;
fig_EKF_xerror_yadata=figure();
plot_errors(tvec,nan(size(xsigmas)),xsigmas, xunits,[40,50,.5,20,20,1])
sgtitle('State Error Estimate, LKF for ydata','FontSize',14, 'Interpreter','latex')
legend('','+-$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')
saveas(fig_EKF_xerror_yadata, 'images/EKF_xerror_ydata_pt6.png')

% LKF
[Pk,dyhat,dxhat,innovation,Sk,xsigmas,ysigmas] = LKF(ydata,dxhat0_LKF,P0_LKF,du0_LKF,tvec,Q_LKF,R_LKF,Abar,Bbar,H,xnom);
yhat_LKF = dyhat + ynom;

% measurements + corresponding KF estimation
fig_LKF_yadata=figure();
plot_KF(tvec, ydata, yhat_LKF, ysigmas, yunits, wrap_indices_y)
sgtitle('Measurements, LKF for ydata','FontSize',14, 'Interpreter','latex')
legend('Ground truth', 'KF output','FontSize',12,'Interpreter','latex')
saveas(fig_LKF_yadata, 'images/LKF_ydata_pt6.png')

% plot innovation vector and 2*sigma_y
fig_LKF_error_yadata=figure();
plot_errors(tvec(2:end),innovation(:,2:end),ysigmas(:,2:end), xunits)
sgtitle('Measurement Error Estimate, LKF for ydata','FontSize',14, 'Interpreter','latex')
saveas(fig_LKF_error_yadata, 'images/LKF_error_ydata_pt6.png')

% plot 2*sigma_x
fig_LKF_xerror_yadata=figure();
plot_errors(tvec,nan(size(xsigmas)),xsigmas, xunits,[40,50,.5,20,20,1])
sgtitle('State Error Estimate, LKF for ydata','FontSize',14, 'Interpreter','latex')
legend('','+-$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')
saveas(fig_LKF_xerror_yadata, 'images/LKF_xerror_ydata_pt6.png')

% state errors for both KFs plotted together
fig_KF_xerror_yadata=figure();
plot_both_states(tvec,xsigmas,xsigmas_EKF, xunits, wrap_indices_x,[40,40,.7,15,15,1])
sgtitle('State Error Estimate for LKF and EKF','FontSize',14, 'Interpreter','latex')
legend('LKF +-$2\sigma$ Error Bound','EKF +-$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')
saveas(fig_KF_xerror_yadata, 'images/all_KF_xerror_ydata_pt6.png')

% Plotting the states themselves
fig_state_est_yadata = figure();

%plot_states(tvec,xhat,xunits,wrap_indices_x)
plot_both_states(tvec,xhat,(dxhat+xnom), xunits, wrap_indices_x)
sgtitle('State Estimation from ydata','FontSize',14, 'Interpreter','latex')
legend('EKF', 'LKF','FontSize',12, 'Interpreter', 'latex');

saveas(fig_state_est_yadata, 'images/state_est.png')

%% Ode45 Function


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

%% CT and DT Matrices Function

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
            (x5-x2)/abv (x1-x4)/abv 0 (x2-x5)/abv (x4-x1)/abv -1; ... % I changed the 6th column to -1 (it used to be 0)
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

%% Plotting Functions

function plot_states(tarr,state,ylabels,wrap_indices,line_type)
% plot the states using subplots
    for iw = 1:length(wrap_indices)
        state(wrap_indices(iw),:)= mod(state(wrap_indices(iw),:)+pi,2*pi)-pi;  
    end
    lw = 2;
    if not(exist('line_type','var'))
        line_type = '-';
    end

    for i=1:size(state,1)
        subplot(size(state,1),1,i)
        hold on
        plot(tarr,state(i,:),line_type,'linewidth',lw)%,'Color','blue','LineWidth',1.5)
        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
        grid on
        hold off
    end

end

% ------ Another Plotting Function ---------
function plot_both_states(tarr,state1,state2,ylabels,wrap_indices,lims)
% plot the states using subplots
    for iw = 1:length(wrap_indices)
        state1(wrap_indices(iw),:)= mod(state1(wrap_indices(iw),:)+pi,2*pi)-pi;  
        state2(wrap_indices(iw),:)= mod(state2(wrap_indices(iw),:)+pi,2*pi)-pi;  
    end

    for i=1:size(state1,1)
        subplot(size(state1,1),1,i)
        plot(tarr,state1(i,:),'Color','blue','LineWidth',1.5)
        hold on
        plot(tarr,state2(i,:),'--','Color','red','LineWidth',1.5)
        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
        grid on
        if exist('lims','var')
            ylim([-lims(i) lims(i)])
        end
        
    end

end

% ------ Error Plotting Function
function plot_errors(tarr,errors,sigmas,ylabels,lims)
% plot the states using subplots
    
    for i=1:size(errors,1)
        subplot(size(errors,1),1,i)
        plot(tarr,errors(i,:),'Color','red','LineWidth',1.5)
        hold on
        plot(tarr,sigmas(i,:),'--','Color','black','LineWidth',1.5)
        hold on
        plot(tarr,-1*sigmas(i,:),'--','Color','black','LineWidth',1.5)
        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
        grid on
        if exist('lims','var')
            ylim([-lims(i) lims(i)])
        end
        
    end
    legend('Error','$2\sigma$ Error Bound','$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')

end

% ------ KALMAN FILTER PLOTTING FUNCTION -------

function plot_KF(tarr,sim_state,KF_state, sigmas, ylabels,wrap_indices,real_data)
if not(exist('real_data','var'))
    truth_label = 'Simulated ground truth';
else 
    truth_label = 'Ground truth';
end

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
        % hold on
        % plot(tarr, KF_state(i,:) + sigmas(i,:),'--','Color','black','LineWidth',1.5)
        % hold on
        % plot(tarr, KF_state(i,:) - sigmas(i,:),'--','Color','black','LineWidth',1.5)

        ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
        xlabel('Time (s)','FontSize',12, 'Interpreter','latex')

        grid on

        legend(truth_label, 'KF output','$2\sigma$ Error Bound','$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')

       % legend('Simulated ground truth', 'KF output','$2\sigma$ Error Bound','$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')

    end

end

function plot_EKF(tarr,sim_state,KF_state,ylabels,wrap_indices,real_data)
if not(exist('real_data','var'))
    truth_label = 'Simulated ground truth';
else 
    truth_label = 'Ground truth';
end

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
        legend(truth_label, 'KF output','FontSize',12, 'Interpreter','latex')
    end

end

function plot_error(tarr,error,sigmas_all, ylabels)
% plot the states using subplots
    % for iw = 1:length(wrap_indices)
    %     sim_state(wrap_indices(iw),:)= mod(sim_state(wrap_indices(iw),:)+pi,2*pi)-pi;  
    %     KF_state(wrap_indices(iw),:)= mod(KF_state(wrap_indices(iw),:)+pi,2*pi)-pi; 
    % end
    if size(sigmas_all,4) ==1
        plot_ind = 1;
    else
        plot_ind = 5;
    end
    

        for i=1:size(error,1)
            subplot(size(error,1),1,i)
            sigma_plus = squeeze(sigmas_all(i,i,:,plot_ind ))';%error(i,:) + squeeze(sigmas_all(i,i,:,plot_ind ))';
            sigma_minus = -1*squeeze(sigmas_all(i,i,:,plot_ind ))';%error(i,:) - squeeze(sigmas_all(i,i,:,plot_ind ))';
            hold on
            plot(tarr,error(i,:),'Color','blue','LineWidth',.5)
            plot(tarr, sigma_plus,"Color",'r', 'LineStyle', '--')
            plot(tarr, sigma_minus, 'Color','r', 'LineStyle', '--')
            hold off
            ylabel(ylabels{i},'FontSize',12, 'Interpreter','latex')
            xlabel('Time (s)','FontSize',12, 'Interpreter','latex')
            grid on
            legend('KF output error','$2\sigma$ Error Bound','$2\sigma$ Error Bound','FontSize',12,'Interpreter','latex')
            sgtitle('State Error Estimate','FontSize',14, 'Interpreter','latex')
        end

end


%% Observations Calculation Functions
% ------- with noise --------
function y = calc_obs_from_state(x,vtilde)
    x1 = x(1,:);
    x2 = x(2,:);
    x3 = x(3,:);
    x4 = x(4,:);
    x5 = x(5,:);
    x6 = x(6,:);
       
    y = [atan2( (x5-x2) , (x4-x1) ) - x3;...
             sqrt( (-x4+x1).^2 + (-x5+x2).^2 );...
             atan2( (x2-x5) , (x1-x4) ) - x6;...
             x4;...
             x5]+vtilde;
end

% ------- without noise --------
function y = compute_Y(x)
    % Extract state variables
    x1 = x(1); % xi_g
    x2 = x(2); % eta_g
    x3 = x(3); % theta_g
    x4 = x(4); % xi_a
    x5 = x(5); % eta_a
    x6 = x(6); % theta_a

    % Compute measurements
     y = [atan2( (x5-x2) , (x4-x1) ) - x3;...
             sqrt( (-x4+x1).^2 + (-x5+x2).^2 );...
             atan2( (x2-x5) , (x1-x4) ) - x6;...
             x4;...
             x5];
end

%% Dynamic Matrices Functions
% -----------DYNAMIC ABAR, BBAR, HBAR MATRICES FOR EKF----------------

function Abar = compute_Abar(x, u)
    % Extract state variables
    x3 = x(3); % theta_g
    x6 = x(6); % theta_a
    
    % Extract input variables
    u1 = u(1); % v_g
    u2 = u(2); % phi_g
    u3 = u(3); % v_a
    
    % Define Abar matrix
    Abar = [0 0 -u1*sin(x3)  0  0  0;...
        0 0  u1*cos(x3)  0  0  0;...
        0 0  0  0  0  0;...
        0 0  0  0  0 -u3*sin(x6);...
        0 0  0  0  0 u3*cos(x6);...
        0 0  0  0  0  0];
end



function Bbar = compute_Bbar(x, u)
    % Extract state variables
    x3 = x(3); % theta_g
    x6 = x(6); % theta_a
    L = .5;
    
    % Extract input variables
    u1 = u(1); % v_g
    u2 = u(2); % phi_g
    
    % Define Bbar matrix
    Bbar = [cos(x3)  0  0  0;...
        sin(x3)  0  0  0;...
        (1/L)*tan(u2) u1/(L*(sec(u2)^2))  0  0;... 
        0  0 cos(x6)  0;...
        0  0 sin(x6)  0;...
        0  0 0  1];
end

function Cbar = compute_Cbar(x)
    % Extract state variables
    x1 = x(1); % xi_g
    x2 = x(2); % eta_g
    x3 = x(3); % theta_g
    x4 = x(4); % xi_a
    x5 = x(5); % eta_a
    x6 = x(6); % theta_a

    % Precompute terms
    abv = (x4 - x1)^2 + (x5 - x2)^2; % Distance squared between ground and air vehicles

    % Define Cbar matrix
    Cbar = [(x5-x2)/abv  (x1-x4)/abv  -1  (x2-x5)/abv  (x4-x1)/abv  0;...
        (x1-x4)/sqrt(abv)  (x2-x5)/sqrt(abv)  0  (x4-x1)/sqrt(abv)  (x5-x2)/sqrt(abv)  0;...
        (x5-x2)/abv  (x1-x4)/abv  0  (x2-x5)/abv  (x4-x1)/abv  -1;...
        0  0  0  1  0  0;...
        0  0  0  0  1  0];

    
end





function [Pk_all, yhat,xhat_p,innovation,Sk_collect,sigmas_collect,xsigmas,ysigmas] = EKF(y_truth_sim,xhat0,Pk0,unom,L,tarr,Q,R)
% Initialize EKF
    dt = tarr(2)-tarr(1);
    xhat_p = zeros(6, length(tarr)); 
    xhat_p(:, 1) = xhat0; % Initial state estimate

    Pk_all = zeros(6, 6, length(tarr)); 
    Pk_all(:, :, 1) = Pk0;

    yhat = zeros(5,length(tarr));
    innovation = zeros(5,length(tarr));
    Sk_collect = zeros(5,5,length(tarr));

    sigmas_collect = zeros(6,6,length(tarr));
    P0 = Pk0;
    sigmas_collect(:,:,1) = diag([2*sqrt(P0(1,1)) 2*sqrt(P0(2,2)) 2*sqrt(P0(3,3)) 2*sqrt(P0(4,4)) 2*sqrt(P0(5,5)) 2*sqrt(P0(6,6))]);
    xsigmas = zeros(6,length(tarr));
    ysigmas = zeros(5,length(tarr));

    for k=2:length(tarr)
        % Calculate Jacobians for each time step using current state
        Abar_k = compute_Abar(xhat_p(:, k-1), unom(:, k-1));
        Bbar_k = compute_Bbar(xhat_p(:, k-1), unom(:, k-1)); 

        Ftild_k = eye(6) + dt * Abar_k; % State transition matrix
        Gtild_k = dt * Bbar_k; % Input matrix.
        Omegatild_k = dt * eye(6);  

        % Compute Q_k dynamically
        Q_k = Q; 
    
        %EKF prediction
        % From lecture notes: assume wk=0
        my_ode = @(t, y) NL_ode(t, y, unom(1, k-1), unom(2, k-1), unom(3, k-1), unom(4, k-1),zeros(3),zeros(3), L); 
        [~, xhat_minus] = ode45(my_ode, [tarr(k-1), tarr(k)], xhat_p(:, k-1));

        xhat_minus = xhat_minus(end, :)'; % Take last output as predicted state
        
    
        % covariance 
        Pk_minus = Ftild_k * Pk_all(:, :, k-1) * Ftild_k' + Omegatild_k * Q_k * Omegatild_k';
   

        % innovation
        predicted_y = compute_Y(xhat_minus);
        Htild_k = compute_Cbar(xhat_minus);

        ey_k = y_truth_sim(:, k) - predicted_y; % Actual measurements - predicted  
        padding = 1e-6; % tiny tolerance
        ey_k(1) = mod(ey_k(1) + pi + padding, 2*pi) - pi - padding;
        ey_k(3) = mod(ey_k(3) + pi + padding, 2*pi) - pi - padding;
        
        % correction
        Sk = Htild_k * Pk_minus * Htild_k' + R;
        Skval = 0.5*(Sk + Sk');

        Kk = Pk_minus * Htild_k' / Sk;

        xhat_plus = xhat_minus + Kk * ey_k;
        Pk_plus = (eye(6) - Kk * Htild_k) * Pk_minus;


        
        ysigma1 = 2*sqrt(Skval(1,1));
        ysigma2 = 2*sqrt(Skval(2,2));
        ysigma3 = 2*sqrt(Skval(3,3));
        ysigma4 = 2*sqrt(Skval(4,4));
        ysigma5 = 2*sqrt(Skval(5,5));

        % extract 2sigma values
        xsigma1 = 2*sqrt(Pk_plus(1,1));
        xsigma2 = 2*sqrt(Pk_plus(2,2));
        xsigma3 = 2*sqrt(Pk_plus(3,3));
        xsigma4 = 2*sqrt(Pk_plus(4,4));
        xsigma5 = 2*sqrt(Pk_plus(5,5));
        xsigma6 = 2*sqrt(Pk_plus(6,6));
        


        % Save results
        xhat_p(:, k) = xhat_plus;
        Pk_all(:, :, k) = Pk_plus;
        yhat(:,k) = predicted_y;
        innovation(:,k) = ey_k;
        Sk_collect(:,:,k) = Skval;
        sigmas_collect(:,:,k) = diag([2*sqrt(Pk_plus(1,1)) 2*sqrt(Pk_plus(2,2)) 2*sqrt(Pk_plus(3,3)) 2*sqrt(Pk_plus(4,4)) 2*sqrt(Pk_plus(5,5)) 2*sqrt(Pk_plus(6,6))]);
        xsigmas(:,k) = [xsigma1;xsigma2;xsigma3;xsigma4;xsigma5;xsigma6];
        ysigmas(:,k) = [ysigma1;ysigma2;ysigma3;ysigma4;ysigma5];

    end
end

function [Pk,dyhat,dxhat,innovation,Sk,xsigmas,ysigmas] = LKF(y_truth_sim,dxhat0_LKF,Pk0_LKF,du0_LKF,tarr,Q_LKF,R,Abar,Bbar,H,xnom)

    
    % ----------------
    % Use the simulated "truth" measurements as inputs to KF

    % initialize KF
    dxhat = zeros(6,length(tarr)); 
    dxhat(:,1) = dxhat0_LKF;

    dyhat = zeros(5,length(tarr));

    du = zeros(4,length(tarr));
    du(:,1) = du0_LKF;

    % initial state covariance matrix
    P0 = Pk0_LKF;
    Pk = zeros(6,6,length(tarr));
    Pk(:,:,1) = P0;

    xsigmas = zeros(6,length(tarr));
    ysigmas = zeros(5,length(tarr));
    xsigmas(:,1) = [2*sqrt(P0(1,1)); 2*sqrt(P0(2,2)); 2*sqrt(P0(3,3)); 2*sqrt(P0(4,4)); 2*sqrt(P0(5,5)); 2*sqrt(P0(6,6))];

    dt = tarr(2)-tarr(1);

    Q_KF = Q_LKF;
    R_KF = R;
    for k = 2:length(tarr)
        % must re-calculate F, G, H, Omega at each timestep! These are not
        % the same as the F, G, H Jacobians calculated in part 1.
        Ftild_k = eye(6) + dt*Abar(:,:,k-1);
        Htild_k = H(:,:,k); % H tilde = H bar = C bar
        Gtild_k = dt*Bbar(:,:,k-1);
        Gamma = eye(6);
        Omegatild_k = dt*Gamma;

        % prediction step
        dxhat_minus = Ftild_k*dxhat(:,k-1) + Gtild_k*du(:,k-1);
        Pk_minus = Ftild_k*Pk(:,:,k-1)*Ftild_k' + Omegatild_k*Q_KF*Omegatild_k';
       
        Skval = Htild_k*Pk_minus*Htild_k' + R_KF;
        Skval = 0.5*(Skval + Skval');

        % gain K
        K = Pk_minus*Htild_k' * inv(Htild_k*Pk_minus*Htild_k' + R_KF);

        % correction step
        Pk_plus = (eye(6) - K*Htild_k)*Pk_minus;
        pretend_y_data = y_truth_sim(:,k); % this is y_k+1
        % y* contains NO NOISE
        predicted_y = compute_Y(xnom(:,k)); % this is y*_k+1
        dy = pretend_y_data-predicted_y;
        innovation_k = dy-Htild_k*dxhat_minus; 
        % angle wrapping - rows 1 and 3 of the innovation vector
            % due to angles being very close to pi or -pi due to atan2
            % we want to shift the angle range to be 0 to 2pi, then
            % subtract pi to eliminate discontinuities near +/-pi
        padding = 1e-6; % tiny tolerance
        innovation_k(1) = mod(innovation_k(1) + pi + padding, 2*pi) - pi - padding;
        innovation_k(3) = mod(innovation_k(3) + pi + padding, 2*pi) - pi - padding;
        dxhat_plus = dxhat_minus + K*(innovation_k);

        % results of the state & covariance
        dxhat(:,k) = dxhat_plus;
        Pk(:,:,k) = Pk_plus;
        dyhat(:,k) = Htild_k*dxhat_minus;
        Sk(:,:,k) = Skval;
        innovation(:,k) = innovation_k;

        % extract 2sigma values
        xsigma1 = 2*sqrt(Pk(1,1,k));
        xsigma2 = 2*sqrt(Pk(2,2,k));
        xsigma3 = 2*sqrt(Pk(3,3,k));
        xsigma4 = 2*sqrt(Pk(4,4,k));
        xsigma5 = 2*sqrt(Pk(5,5,k));
        xsigma6 = 2*sqrt(Pk(6,6,k));
        xsigmas(:,k) = [xsigma1;xsigma2;xsigma3;xsigma4;xsigma5;xsigma6];

        ysigma1 = 2*sqrt(Sk(1,1,k));
        ysigma2 = 2*sqrt(Sk(2,2,k));
        ysigma3 = 2*sqrt(Sk(3,3,k));
        ysigma4 = 2*sqrt(Sk(4,4,k));
        ysigma5 = 2*sqrt(Sk(5,5,k));
        ysigmas(:,k) = [ysigma1;ysigma2;ysigma3;ysigma4;ysigma5];

    end
end

