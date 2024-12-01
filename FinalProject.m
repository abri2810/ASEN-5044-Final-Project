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

% Observability
%ob = [H; H*F; H*F^2; H*F^3; H*F^4; H*F^5];
%rank(ob);
% The observability matrix has full column rank, so it is observable

% Stability
%eig(F);
% The system is marginally stable because all eigenvalues lie on the unit
% circle.

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

% initialize stuff 
P0 = eye(6); % initial state covariance matrix (IDK WHAT TO PUT HERE SO I MADE IT IDENTITY)
MC_num = 100; % number of monte carlo simulations
NEES = zeros(MC_num,length(tarr));
NIS = zeros(MC_num,length(tarr));
dxhat_all = zeros(6,length(tarr),MC_num);
noisy_y = zeros(5, length(tarr),MC_num);


% ----- super iffy on the code below ----------

for m = 1:MC_num % for each MC iteration
    % simulate noisy measurement for each iteration
        % create random measurement noise
        v = mvnrnd([0;0;0;0;0],R,length(tarr));
        % simulate measurement 
        noisy_y(:,:,m) = ydata + v' ;
            % this gives a new noisy measurement for each monte carlo
            % iteration. this will be used as input to the KF to produce
            % state estimates and predicted measurements

    % initialize KF
    dxhat0 = deltx0; % initial PERTURBATION state estimate
    P = eye(6); % initial state covariance matrix (IDK WHAT TO PUT HERE SO I MADE IT IDENTITY)
    dxhat = zeros(6,length(tarr));
    dxhat(:,1) = xhat0;
    du0 = [0;0;0;0;0;0]; % ADD REAL NUMBERS TO THIS LATER
    du = zeros(6,length(tarr));

    for k = 2:length(tarr) % for each timestep k 
        F_t = F(:,:,k-1); % time-varying F matrix
        H_t = H(:,:,k-1); % time-varying H matrix
        G_t = G(:,:,k-1); % time varying G matrix

        % prediction step
        dxhat_minus = F_t*dxhat(:,k-1) + G_t*du(:,k-1);

    end
    dxhat_all(:,:,m) = dxhat;
end


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
    Abar = [0 0 -u1*sin(x3) 0 0 0; ...
            0 0 u1*cos(x3) 0 0 0;...
            0 0 0 0 0 0; ...
            0 0 0 0 0 -u3*sin(x6); ...
            0 0 0 0 0 u3*cos(x6); ...
            0 0 0 0 0 0];
    
    
    Bbar = [cos(x3) 0 0 0; ...
            sin(x3) 0 0 0; ...
            (1/L)*tan(u2) u1/L*(sec(u2))^2 0 0; ...
            0 0 cos(x6) 0; ...
            0 0 sin(x6) 0; ...
            0 0 0 1];
    
    abv = (x4-x1)^2 + (x5-x2)^2;
    
    Cbar = [(x5-x2)/abv (x1-x4)/abv -1 (x2-x5)/abv (x4-x1)/abv 0; ...
            (x1-x4)/sqrt(abv) (x2-x5)/sqrt(abv) 0 (x4-x1)/sqrt(abv) (x5-x2)/sqrt(abv) 0; ...
            (x5-x2)/abv (x1-x4)/abv 0 (x2-x5)/abv (x4-x1)/abv 0; ...
            0 0 0 1 0 0; ...
            0 0 0 0 1 0];


    Dbar = zeros(5,4);
    
    
    % Part I, Problem 2.
    
    z = [Abar Bbar; zeros(4,6) zeros(4)];
    ez = expm(z*dt);
    
    F(:,:,ti) = ez(1:6, 1:6);
    G(:,:,ti) = ez(1:6, 7:10);
    H(:,:,ti) = Cbar;
    M(:,:,ti) = Dbar;
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