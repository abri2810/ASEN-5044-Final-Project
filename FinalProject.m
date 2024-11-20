% 5044 Final Project
% Sarah Luettgen, Abby Rindfuss, and Lisa Ventura
% Cooperative Location

% Housekeeping
clear; 

% Part I
% Problem 1.

L = .5;
%phi_g is between -5*pi/12 to 5*pi/12
%v_gmax = 3;
%omega_g is between -pi/6 to pi/6
%v_a is between 10 and 20
xi_g0 = 10;
eta_g0 = 0;
theta_g = pi/2;
v_g = 2;
phi_g = -pi/18;

xi_a0 = -60;
eta_a0 = 0;
theta_a = -pi/2;
v_a = 12;
omega_a = pi/25;

dt = .1;

u1 = [v_g; phi_g];
u2 = [v_a; omega_a];
x1 = xi_g0;
x2 = eta_g0;
x3 = theta_g;
x4 = xi_a0;
x5 = eta_a0;
x6 = theta_a;

Abar = [0 0 -u1(1,:)*sin(x3) 0 0 0; 0 0 u1(2,:)*cos(x3) 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 -u2(1,:)*sin(x6); ...
        0 0 0 0 0 u2(2,:)*cos(x6); 0 0 0 0 0 0];

Bbar = [cos(x3) 0 0 0; sin(x3) 0 0 0; (1/L)*tan(u2(1,:)) (u1(1,:)/L)*sec(u2(2,:)).^2 0 0; ...
        0 0 cos(x6) 0; 0 0 sin(x6) 0; 0 0 0 1];

abv = (x4-x1)^2 + (x5-x2)^2;
Cbar = [(x5-x2)/abv (x1-x4)/abv -1 (x2-x5)/abv (x4-x1)/abv 0; ...
        (x1-x4)/sqrt(abv) (x2-x5)/sqrt(abv) 0 (x4-x1)/sqrt(abv) (x5-x2)/sqrt(abv) 0; ...
        (x5-x2)/abv (x1-x4)/abv 0 (x2-x5)/abv (x4-x1)/abv 0; ...
        0 0 0 1 0 0; ...
        0 0 0 0 1 0];

Dbar = zeros(5,4);


%% Part 2

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

