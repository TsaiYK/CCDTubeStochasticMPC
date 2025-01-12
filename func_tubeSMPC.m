function [V_N,Size_Sigma_e] = func_tubeSMPC(xp,K,SystemType)

%% Define system and parameters
System = LinearSystemDef;
switch SystemType
    case '2D_numerical'
        System.numerical_example(xp,K,1);
        mean_epsilon = 0; std_epsilon = 0.1;
        % where epsilon_k follows a normal distribution with 0 mean and std_epsilon std
    case '4D_satellite'
        System.satellite_dynamics(xp,K,1);
        mean_epsilon = 0; std_epsilon = 0.25;
        % where epsilon_k follows a normal distribution with 0 mean and std_epsilon std
end
Ak_inf = System.A_K^1e3;

%% Disturbance: x_{k+1} = A*x_k + B*u_k + w_k, where w_k = B*epsilon_k

mu_w = System.B*mean_epsilon;
Sigma_w = diag(System.B*std_epsilon').^2;

Constr = ConstraintTightening;
Constr.nominal_constr_para(System,K,mu_w,Sigma_w);

mpc = ModelPredictiveControl(System.mysys, Constr.Xc, Constr.Uc,...
    Constr.Xc_bar, Constr.Uc_bar, System.N_horizon);

% if isempty(mpc.XMPI)
%     V_N = NaN;
%     Size_Sigma_e = NaN;
% else
x_curr = System.iniCon;
infeasible_mpc = false;
for i = 1:System.N
    %% solve the MPC problem
    u_next = mpc.solve(x_curr);
    if isempty(u_next)
        infeasible_mpc = true;
        break
    else
        %% collect the nominal control inputs and states
        x_nom = mpc.solution_cache.x_nominal_seq;
        u_nom = mpc.solution_cache.u_nominal_seq;
        
        %% Compute loss function
        L(i) = x_nom(:,1)'*System.Q*x_nom(:,1)+u_nom(:,1)'*System.R*u_nom(:,1);
        if i==System.N % include terminal cost
            L(i+1) = x_nom(:,2)'*System.mysys.P*x_nom(:,2);
        end
        
        % update the current state (predicted nominal state will be the current state at the next time)
        x_curr = x_nom(:,2);
    end
end
if infeasible_mpc
    V_N = NaN;
    Size_Sigma_e = NaN;
elseif sum(find(isinf(Ak_inf)),'all') || sum(find(isnan(Ak_inf)),'all')
    V_N = NaN;
    Size_Sigma_e = NaN;
elseif norm(eig(Ak_inf))>1e-2
    V_N = NaN;
    Size_Sigma_e = NaN;
    %     elseif ~(x_curr<=mpc.XMPI) % check if the final state is in the MPI set
    %         V_N = NaN;
    %         Size_Sigma_e = NaN;
else
    clear('mpc')
    % Computing the norm of the eigenvalues of the Sigma_e_{\infty}
    Sigma_e_N = Constr.Sigma_e{1};
    for i = 2:1e4
        Sigma_e_N = Sigma_e_N+System.A_K^(i-1)*Constr.Sigma_e{1}*(System.A_K^(i-1))';
        if norm(eig(System.A_K^(i-1)*Constr.Sigma_e{1}*(System.A_K^(i-1))'))<1e-10
            break
        end
    end
    if i<=1e4
        V_N = sum(L);
        Size_Sigma_e = norm(eig(Sigma_e_N));
    else
        % if Sigma_e_N cannot converge with N<=1e4, it cannot be counted
        % as a feasible design.
        V_N = NaN;
        Size_Sigma_e = NaN;
    end
    if abs(Size_Sigma_e)>0.5
        V_N = NaN;
        Size_Sigma_e = NaN;
    end
end
% end
