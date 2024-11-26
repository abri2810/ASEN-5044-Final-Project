% 5044 Final Project
% Sarah Luettgen, Abby Rindfuss, and Lisa Ventura
% Cooperative Location

% Housekeeping
clear; 
clc
close all

%% Part I, Problem 1.

% ---- nominal conditions -----

L = 0.5; % UGV base length
    % phi_g is between -5*pi/12 to 5*pi/12
    % v_gmax = 3;
    % omega_g is between -pi/6 to pi/6
    % v_a is between 10 and 20
xi_g0 = 10;
eta_g0 = 0;
theta_g = pi/2;
v_g = 2; % UGV nominal speed (m/s)
phi_g = -pi/18; % UGV nominal steering angle (rad)

xi_a0 = -60;
eta_a0 = 0;
theta_a = -pi/2;
v_a = 12; % UAV nominal s peed (m/s)
omega_a = pi/25; % UAV nominal turning rate (rad/s)

% delta T sampling rate
dt = .1;

% nominal u(t) vector components 
u1 = v_g;
u2 = phi_g;
u3 = v_a;
u4 = omega_a;

% nominal x(t) vector components 
x1 = xi_g0;
x2 = eta_g0;
x3 = theta_g;
x4 = xi_a0;
x5 = eta_a0;
x6 = theta_a;


% ------- Jacobians -------

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


%% Part I, Problem 2.

z = [Abar Bbar; zeros(4,6) zeros(4)];
ez = expm(z*dt);

F = ez(1:6, 1:6);
G = ez(1:6, 7:10);
H = Cbar;
M = Dbar;


% Observability
ob = [H; H*F; H*F^2; H*F^3; H*F^4; H*F^5];
rank(ob);
% The observability matrix has full column rank, so it is observable

% Stability
eig(F);
% The system is marginally stable because all eigenvalues lie on the unit
% circle.

%% Part I, Problem 3.

tf = 100; %seconds

% initialize state vectors
tarr = 0:dt:(tf); % t vector
du = zeros(4,length(tarr)); % du vector
deltx0 = [0; 1; 0; 0; 0; 0.1]; %dx0
% deltx0 = [x1; x2; x3; x4; x5; x6];
d_state = nan(6,length(tarr)); %dx
d_state(:,1) = deltx0;

% nominal solution
vg_nom = 2; %m/s
phi_nom = -pi/18; %rad
va_nom = 12; %m/s
omegaa_nom = pi/25; %rad/s
thetag_dot_nom = v_g/L*tan(phi_nom);
theta_g_nom = theta_g+thetag_dot_nom*tarr;
thetaa_dot_nom = omegaa_nom;
theta_a_nom = theta_a+thetaa_dot_nom*tarr;

xnom = [xi_g0+v_g/thetag_dot_nom*sin(theta_g_nom)-v_g/thetag_dot_nom*sin(theta_g);...
        eta_g0-v_g/thetag_dot_nom*cos(theta_g_nom)-v_g/thetag_dot_nom*cos(theta_g);...
        theta_g_nom;...
        xi_a0+v_a/thetaa_dot_nom*sin(theta_a_nom)-v_a/thetaa_dot_nom*sin(theta_a);...
        eta_a0-v_a/thetaa_dot_nom*cos(theta_a_nom)-v_a/thetaa_dot_nom*cos(theta_a);...
        theta_a_nom];

% solve for dx_k+1
for i = 2:length(tarr)
    d_state(:,i) = F*d_state(:,i-1) + G*du(:,i-1);
end

full_state = xnom + d_state; % x = x_nom + dx

% full state solved with ODE45
my_ode = @(t,y) NL_ode(t,y,vg_nom,phi_nom,va_nom,omegaa_nom,[0;0;0],[0;0;0],L);
[t,yarr] = ode45(my_ode,tarr,[x1;x2;x3;x4;x5;x6]+deltx0);
yarr=yarr';
yarr(3,:) = mod(yarr(3,:)+pi,2*pi)-pi;
yarr(6,:) = mod(yarr(6,:)+pi,2*pi)-pi;

% plot perturbations
figure
for i=1:size(full_state,1)
    subplot(size(full_state,1),1,i)
    %plot(tarr,d_state(i,:),DisplayName="linear")
    plot(tarr,d_state(i,:))
end

% full state
figure
for i=1:size(full_state,1)
    subplot(size(full_state,1),1,i)
    %plot(tarr,d_state(i,:),DisplayName="linear")
    plot(tarr,full_state(i,:))
end


%% Part II, Problem 4. 

coopData = load('cooplocalization_finalproj_KFdata.mat');
Q = coopData.Qtrue;
R = coopData.Rtrue;
ydata = coopData.ydata;


%% Functions
function yd = NL_ode(t,y,vg,phi,va,wa,w_tild_g,w_tild_a,L)
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
