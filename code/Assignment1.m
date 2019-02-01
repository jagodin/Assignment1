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
dim_x = 200e-9;
dim_y = 100e-9;

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
% TODO: optimize calculations

num_e = 10; % number of electrons

% initialize x and y position of electrons
[x_vec, y_vec] = initPosition(num_e, dim_x, dim_y); 

% initialize x and y velocity of electrons
[vx_vec, vy_vec] = initVelocity(num_e, Vth);

% initialize time variables
t = 0;
steps = 100;
t_step = max(dim_x, dim_y)/(50*Vth);
t_final = steps*t_step;
t_vec = zeros(1,steps+1);

% initialize time vector with time step
t_vec = zeros(1,length(t_vec));
for i=1:length(t_vec)
    t_vec(i) = (i-1)*t_step;
end

colour = hsv(num_e);
j=0;
Temp = zeros(1,length(t_vec));
while t < t_final
    j=j+1;
    % Calculate temp
    Temp(j) = (mean(vx_vec.^2 + vy_vec.^2)*mn)/(2*k);
    
    % Save previous positions
    x_vec_prev = x_vec;
    y_vec_prev = y_vec;
    
    % Calculate new position
    x_vec = x_vec + vx_vec*t_step;
    y_vec = y_vec + vy_vec*t_step;
    
    % Boundary conditions
    for i=1:num_e
        if x_vec(i) < 0 % left boundary, periodic
            x_vec(i) = x_vec(i)+dim_x;
            x_vec_prev(i) = dim_x;
        end
        if x_vec(i) > dim_x % right boundary, periodic
            x_vec(i) = x_vec(i)-dim_x;
            x_vec_prev(i) = 0;
        end
        if y_vec(i) > dim_y % top boundary, reflect
            vy_vec(i) = -vy_vec(i);
            y_vec(i) = 2*dim_y - y_vec(i);
        end
        if y_vec(i) < 0 % bottom boundary, reflect
            vy_vec(i) = -vy_vec(i);
            y_vec(i) = abs(y_vec(i));
        end
    end
    
    % Plot trajectories
    xlabel('x (m)')
    ylabel('y (m)')
    title('Particle Trajectories')
    xlim([0 dim_x])
    ylim([0 dim_y])
    pause(0.1)
    
    for i=1:num_e
        plot([x_vec_prev(i);x_vec(i)],[y_vec_prev(i);y_vec(i)],'color',colour(i,:))
        hold on
    end
    
    t=t+t_step
end

% Plot temp
figure;
plot(t_vec, Temp)
xlabel('time (s)')
ylabel('temp (K)')
ylim([0 500])
xlim([0 t_final])
title('Temperature (K) vs. Time (s)')
pause(0.1)

% ------------------------------------------------------------------------
%
% Part 2: Collisions with Mean Free Path (MFP)
%
% Question 1: Assign a random velocity to each of the particles at the
% start. Use a Maxwell-Boltzmann distribution for each velocity component.
% The average of all speeds will be Vth. Plot distribution in a histogram.
%
% ------------------------------------------------------------------------





% Question 2: Model scattering of electrons using exponential scattering
% probability: P = 1 - exp(-dt/Tmn) where dt is the time since last time
% step. Uf P > rand() then the particle scatters. When the electron
% scatters, re-thermalize its velocities and assign new velocities Vx and
% Vy from the Maxwell-Boltzmann distributions.






