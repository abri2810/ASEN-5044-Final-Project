% 5044 Final Project
% Sarah Luettgen, Abby Rindfuss, and Lisa Ventura
% Cooperative Location

% Housekeeping
clear; 

% Part I
% Problem 1.

syms u1, syms u2, syms u3, syms x1, syms x2, syms x3, syms x4, syms x5, syms x6, syms L;
dt = .1;

Abar = [0 0 -u1*sin(x3) 0 0 0; 0 0 u1*cos(x3) 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 -u3*sin(x6); ...
        0 0 0 0 0 u3*cos(x6); 0 0 0 0 0 0];

Bbar = [cos(x3) 0 0 0; sin(x3) 0 0 0; (1/L)*tan(u2) (u1/L)*sec(u2).^2 0 0; ...
        0 0 cos(x6) 0; 0 0 sin(x6) 0; 0 0 0 1];

abv = (x4-x1)^2 + (x5-x2)^2;
Cbar = [(x5-x2)/abv (x1-x4)/abv -1 (x2-x5)/abv (x4-x1)/abv 0; ...
        (x1-x4)/sqrt(abv) (x2-x5)/sqrt(abv) 0 (x4-x1)/sqrt(abv) (x5-x2)/sqrt(abv) 0; ...
        (x5-x2)/abv (x1-x4)/abv 0 (x2-x5)/abv (x4-x1)/abv 0; ...
        0 0 0 1 0 0; ...
        0 0 0 0 1 0];

Dbar = zeros(5,4);
