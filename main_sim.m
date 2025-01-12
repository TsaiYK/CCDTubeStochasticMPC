clear
clc
close all

addpath(genpath(pwd))
addpath(genpath('src'))
addpath(genpath('tbxmanager'))
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1.5);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',20)

%% Define system and parameters
System = LinearSystemDef;

% 2-D numerical example
xp = 1;
K = [-0.6167,-1.2703];
control_limit = 1;
System.numerical_example(xp, K, control_limit);

% % 4-D satellite system
% xp = cosd(39.14);
% K = [0.182202424222581,-0.438220647924188,-1.109039767665007,-0.978055308689462];
% control_limit = 1;
% System.satellite_dynamics(xp,K,control_limit);

%% Disturbance: x_{k+1} = A*x_k + B*u_k + w_k, where w_k = B*epsilon_k
% where epsilon_k follows a normal distribution with 0 mean and std_epsilon std

% 2-D numerical example
mean_epsilon = 0; std_epsilon = 0.1;

% % 4-D satellite system
% mean_epsilon = 0; std_epsilon = 0.25;

mu_w = System.B*mean_epsilon;
Sigma_w = diag(System.B*std_epsilon').^2;

Constraint = ConstraintTightening;
Constraint.nominal_constr_para(System,K,mu_w,Sigma_w);

mpc = ModelPredictiveControl(System.mysys, Constraint.Xc, Constraint.Uc,...
    Constraint.Xc_bar, Constraint.Uc_bar, System.N_horizon);

%% Dynamic Iteration
savedir_name = './results_mpc_test/'; % saving direction (folder name)
mkdir(savedir_name) % make new folder
n = 1; % number of random samples (if you want to generate the GIFs, just set 1)

x_nom = System.iniCon;
u = cell(1,n);
x = cell(1,n);
for k = 1:n
    x{k}(:,1) = x_nom;
end

w = cell(n,System.N);
for i = 1:n
    for j = 1:System.N
        w{i,j} = [];
        for k = 1:System.nx
            w{i,j} = [w{i,j}; normrnd(mu_w(k),sqrt(Sigma_w(k,k)),[1,System.N_horizon])];
        end
    end
end

% Initialize the state (assuming that there is no disturbance on the initial state)
x_curr = System.iniCon; x_exp = System.iniCon;
x_for_plot = []; u_for_plot = [];
file_name_state = cell(1,System.N);
file_name_ctrl = cell(1,System.N);
file_name_plane = cell(1,System.N);
for i = 1:System.N
    % solve the MPC problem
    mpc.solve(x_curr);
    
    % nominal control and state
    x_nom = mpc.solution_cache.x_nominal_seq;
    u_nom = mpc.solution_cache.u_nominal_seq;
    
    % this loop is to generate states and control inputs with disturbances
    for k = 1:n
        x{k} = x_nom + [zeros(System.nx,1),w{k,i}];
        u{k} = u_nom + K*(x{k}(:,1:end-1)-x_nom(:,1:end-1));
    end
    u_next = u_nom(:,1) + K*(x_curr-x_exp);
    
    % Vectors for plots
    x_for_plot = [x_for_plot,x_curr];
    u_for_plot = [u_for_plot,u_next];
    [prob_vio_state,prob_vio_ctrl] = PlotGenerator(System,Constraint,...
        mpc,x_nom,u_nom,x,u,x_for_plot,u_for_plot,i,n,savedir_name);
    
    % Compute loss function
    L(i) = x_nom(:,1)'*System.Q*x_nom(:,1)+u_nom(:,1)'*System.R*u_nom(:,1);
    if i==System.N
        L(i+1) = x_nom(:,2)'*System.mysys.P*x_nom(:,2);
    end

    % update the current state (predicted nominal state will be the current state at the next time)
    % Use this if you are going to generate the GIF
    x_curr = x{1}(:,2);
    x_exp = x_nom(:,2);
    
    % Use this if you are expecting ideal situation (only nominal states are taken)
    % x_curr = x_nom(:,2);
    % x_exp = x_nom(:,2);
    
    file_name_state{i} = strcat('tmpc_state_seq',number2string(i),'.png');
    file_name_ctrl{i} = strcat('tmpc_ctrl_seq',number2string(i),'.png');
    file_name_plane{i} = strcat('tmpc_plane_seq',number2string(i),'.png');

end

%% Convert the images to GIFs
lps = 'Forever'; delay = 1; delay1 = 1; delay2 = 1;
animated_gif_creator(file_name_state, savedir_name, 'GIF_state.gif', savedir_name, lps, delay, delay1, delay2);
animated_gif_creator(file_name_ctrl, savedir_name, 'GIF_ctrl.gif', savedir_name, lps, delay, delay1, delay2);
animated_gif_creator(file_name_plane, savedir_name, 'GIF_plane.gif', savedir_name, lps, delay, delay1, delay2);

%% Results (Objective values)
V_N = sum(L);
Sigma_e_N = Constraint.Sigma_e{1};
for i = 2:1e4
    Sigma_e_N = Sigma_e_N+System.A_K^(i-1)*Constraint.Sigma_e{1}*(System.A_K^(i-1))';
    if norm(eig(System.A_K^(i-1)*Constraint.Sigma_e{1}*(System.A_K^(i-1))'))<1e-10
        break
    end
end
Size_Sigma_e = norm(eig(Sigma_e_N));
