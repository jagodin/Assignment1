% ELEC 4700
% Assignment 1
%
% Monte-Carlo Modeling of Electron Transport
% 
% Jacob Godin - 100969991
%
% --------------------------------------------
%
% Modeling of the carriers as a population 
% of electrons in an N-type Si semiconductor 
% crystal.
%
% Effective mass of electrons mn = 0.26m0
% Nominal size of the region is 200nm X 100 nm
%
% --------------------------------------------
clear all; clc;

m0 = 9.11e-31;
mn = 0.26*m0;
dim_x = 200;%e-9;
dim_y = 100;%e-9;

k = 1.38064852e-23;

% Part 1: Electron Modelling
%
% Question 1: What is the thermal velocity Vth? Assume T = 300K
%
% Vth = sqrt(vx^2 + vy^2) = sqrt(2kT/mn)

T = 300;
Vth = sqrt(2*k*T/mn);

% Question 2: If the mean time between collisions is Tmn = 0.2ps what is
% the mean free path?
%
% Mean free path = Tmn * Vth

Tmn = 0.2e-12;
Mfp = Tmn * Vth;

% Question 3: Write a pogram that will model the random motion of electrons

num_e = 10; % number of electrons

% initialize x and y position of electrons
[x_vec, y_vec] = initPosition(num_e, dim_x, dim_y); 

% initialize x and y velocity of electrons
[vx_vec, vy_vec] = initVelocity(num_e, Vth);

% initialize time variables
t = 0;
steps = 1000;
t_step = max(dim_x, dim_y)/Vth;
t_final = steps*t_step;
t_vec = zeros(1,steps+1);
for i=1:length(t_vec)
    t_vec(i) = (i-1)*t_step;
end

while t < t_final
    t=t+t_step;
    
    % Boundary conditions
    for i=1:num_e
        if x_vec(i) < 0 % left boundary
            x_vec(i)=x_vec(i)+dim_x;
        end
        if x_vec(i) > dim_x % right boundary
            x_vec(i)=x_vec(i)-dim_x;
        end
        if y_vec(i) > dim_y % top boundary
            
        end
        if y_vec(i) < 0 % bottom boundary
            
        end
    end
    
    % Calculate temp
    
    % Plot temp
    
    % Plot trajectories
    
end