
clear 
close all
clc

%% ALINEA

% simulation initial time step
k_init = 1;
% simulation final time step
k_final = 100;
% entire simulation time length
Tsim = k_final - k_initial + 1;

% disturbance vector: given for the entire simulation period extracted from
% Beats (here I give it a random value)
dist = ones (Tsim,1);

% desired state/occupancy for the link: may be determined from the link's
% capacity
x_des = 30;

% state/occupancy: this is the link's current state that is received from Beats
% (here I give it a random value)
x=10;

% ALINEA's regulator parameter: may be determined by optimizing the
% control input along the prediction horizon
ka = 10;

% initial control input for ALINEA
u_alinea=5;


% Implementing the ALINEA controller within the simulation period 
for i = k_init : k_final
    
    u_alinea = alinea (u_alinea, ka, x_des, x);
    c = 1;
    u_alinea_seq(c) = u_alinea; 
    c = c + 1;
    
    % linking Alinea and a model of the network to create a control input
    % sequence across the prediction horizon
    xp=x;

    for j = i : i + Np - 1
        
        % alinea_model can be any mathematical model of the traffic network
        % or in our case Beats
        xp = alinea_model (u_alinea, xp, dist(j));
        u_alinea = alinea (u_alinea, ka, x_des, xp);
        u_alinea_seq(c) = u_alinea;
        c = c + 1;
        
    end
    
    
    % saving the control inputs produced by ALINEA
    Alinea(i,:) = u_alinea_seq;
    
end

%% Conventional MPC controller

% prediction horizon of the MPC controller (for the entire integrated
% base-parallel control architecture this horizon should be the same)
Np = 5;

% control horizon : to reduce the computation burden (if you don't want it
% set Nc=Np)
Nc=3;

% optimization starting vector: this is received from a base controller
% (here I give it a random value)
u_CMPC = ones (Np, 1);

% Implementing the conventional MPC controller within the simulation period 
for i = k_init : k_final
    
    [u_CMPC, ~, ~, ~] = CMPC (Np, Nc, x, u_CMPC, dist(i:i+Np-1));
    
    % saving the control inputs produced by the conventione MPC controller
    CMPC(i,:) = u_CMPC;
    
end


%% Parameterized MPC controller

% optimization starting vector: this is received from a base controller
% (here I give it a random value)
p_PMPC = ones (Np, 1);

% Implementing the conventional MPC controller within the simulation period 
for i = k_init : k_final
    
    [p_PMPC, ~, ~, ~] = PMPC (Np, Nc, x, p_CMPC, dist(i:i+Np-1));
    
    % saving the control inputs produced by the conventione MPC controller
    PMPC(i,:) = p_CMPC;
    
end
