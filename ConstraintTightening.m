classdef ConstraintTightening < handle
    properties (SetAccess = private)
        eta
        zeta
        mu_e
        Sigma_e
        Sigma_u
        Xc
        Uc
        Xc_bar
        Uc_bar
        Z
        Xf_bar
    end
    
    %% Public functions
    methods (Access = public)
        function obj = nominal_constr_para(obj, sys, K, mu_w, Sigma_w)
            % nominal constraint parameters: eta and zeta such that
            % H*x_nom<=eta; G*u_nom<=zeta
            % mean and covariance matrix of error state: mu_e and Sigma_e
            [obj.eta,obj.zeta,obj.mu_e,obj.Sigma_e] = ...
                obj.find_nom_constr(sys, K, mu_w, Sigma_w);
            
            % covariance matrix of input: Sigma_u
            obj.Sigma_u = zeros(1,length(obj.Sigma_e));
            for i = 1:length(obj.Sigma_e)
                obj.Sigma_u(i) = K*obj.Sigma_e{i}*K';
            end
            
            % original and time-varying constraint sets
            [obj.Xc, obj.Uc, obj.Xc_bar, obj.Uc_bar] = obj.timevarying_constraint_set(sys);
            
            % compute terminal constraint set Xf
            Xf = sys.mysys.compute_MPIset(obj.Xc,obj.Uc);
            
            % obtain the nominal terminal constraint set \bar{X}_f
            obj.Xf_bar = obj.nominal_TerminalSet(Xf);
        end
    end
    
    %% Private functions
    methods (Access = private)
        function [eta,zeta,mu_e,Sigma_e] = find_nom_constr(obj, sys, K, mu_w, Sigma_w)
            %% Define mean and covariance matrix for disturbance w using mu_w and Sigma_w
            mu_e = zeros(sys.nx,sys.N+1);
            mu_e(:,1) = mu_w;
            Sigma_e = cell(1,sys.N+1);
            Sigma_e{1} = Sigma_w;

            for i = 2:(sys.N+1)
                mu_e(:,i) = sys.A_K*mu_e(:,i-1)+mu_w;
                Sigma_e{i} = sys.A_K*Sigma_e{i-1}*sys.A_K'+Sigma_w;
            end

            %% Chance constraints for state and input
            % State constraints
            % Define a new random variable: X_{xi,k}:=H*e_k-h
            mu_X_xi = zeros(size(sys.H,1),sys.N+1);
            std_X_xi = zeros(size(sys.H,1),sys.N+1);

            for k = 1:(sys.N+1)
                for j = 1:size(sys.H,1) % number of state constraints
                    mu_X_xi(j,k) = sys.H(j,:)*mu_e(:,k)-sys.h(j,1);
                    std_X_xi(j,k) = sqrt(sys.H(j,:)*Sigma_e{k}*sys.H(j,:)');
                end
            end

            tol = 0.05; % probability of satisfaction is 1-epsilon
            z_x = norminv(1-tol);
            eta = -(mu_X_xi + z_x*std_X_xi);

            % Input constraints
            % Define a new random variable: X_{u,k}:=G*K*e_k-g
            mu_X_u = zeros(size(sys.G,1),sys.N);
            std_X_u = zeros(size(sys.G,1),sys.N);

            for k = 1:sys.N
                for j = 1:size(sys.G,1) % number of state constraints
                    mu_X_u(j,k) = sys.G(j,:)*K*mu_e(:,k)-sys.g(j,1);
                    std_X_u(j,k) = sqrt(sys.G(j,:)*K*Sigma_e{k}*K'*sys.G(j,:)');
                end
            end

            tol_u = 0.05; % probability of satisfaction is 1-delta
            z_u = norminv(1-tol_u);
            zeta = -(mu_X_u + z_u*std_X_u);
            
        end
        
        function [Xc, Uc, Xc_bar, Uc_bar] = timevarying_constraint_set(obj,sys)
            % Construct the time-varying constraint sets for state and control
            % Xc: original constraint set for state
            % Uc: original constraint set for input
            % Xc_bar: time-varying constraint set for nominal state
            % Uc_bar: time-varying constraint set for nominal input
            xi_max = zeros(sys.nx,sys.N_horizon+1);
            xi_min = zeros(sys.nx,sys.N_horizon+1);
            ui_max = zeros(sys.nu,sys.N_horizon);
            ui_min = zeros(sys.nu,sys.N_horizon);
            xi_max(:,1) = sys.xmax';
            xi_min(:,1) = sys.xmin';
            Xc_bar = cell(1,sys.N_horizon+1);
            Uc_bar = cell(1,sys.N_horizon);
            if sys.nx==2
                Xc_vertex = [xi_max(1,1), xi_max(2,1);...
                    xi_max(1,1), xi_min(2,1);...
                    xi_min(1,1), xi_min(2,1);...
                    xi_min(1,1), xi_max(2,1)];
            elseif sys.nx==4
                Xc_vertex = [xi_max(1,1), xi_max(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_max(1,1), xi_min(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_min(1,1), xi_min(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_min(1,1), xi_max(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_max(1,1), xi_max(2,1), xi_max(3,1), xi_min(4,1);...
                    xi_max(1,1), xi_min(2,1), xi_max(3,1), xi_min(4,1);...
                    xi_min(1,1), xi_min(2,1), xi_max(3,1), xi_min(4,1);...
                    xi_min(1,1), xi_max(2,1), xi_max(3,1), xi_min(4,1);...
                    xi_max(1,1), xi_max(2,1), xi_min(3,1), xi_min(4,1);...
                    xi_max(1,1), xi_min(2,1), xi_min(3,1), xi_min(4,1);...
                    xi_min(1,1), xi_min(2,1), xi_min(3,1), xi_min(4,1);...
                    xi_min(1,1), xi_max(2,1), xi_min(3,1), xi_min(4,1);...
                    xi_max(1,1), xi_max(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_max(1,1), xi_min(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_min(1,1), xi_min(2,1), xi_max(3,1), xi_max(4,1);...
                    xi_min(1,1), xi_max(2,1), xi_max(3,1), xi_max(4,1)];
            else
                error('Sorry! Now we only accept the systems with 2 or 4 dimensions. Please make sure sys.nx = 2 or 4')
            end
            Xc = Polyhedron(Xc_vertex);
            Uc = Polyhedron([sys.umax; sys.umin]);
            Xc_bar{1} = Xc;
            for k = 1:sys.N_horizon
                xi_max(:,k+1) = obj.eta(1:2:end,k);
                xi_min(:,k+1) = -obj.eta(2:2:end,k);
                ui_max(:,k) = obj.zeta(1,k);
                ui_min(:,k) = -obj.zeta(2,k);

                % constraints on state Xc and input Uc
                if sys.nx==2
                    Xc_vertex = [xi_max(1,k+1), xi_max(2,k+1);...
                        xi_max(1,k+1), xi_min(2,k+1);...
                        xi_min(1,k+1), xi_min(2,k+1);...
                        xi_min(1,k+1), xi_max(2,k+1)];
                elseif sys.nx==4
                    Xc_vertex = [xi_max(1,k+1), xi_max(2,k+1), xi_max(3,k+1), xi_max(4,k+1);...
                        xi_max(1,k+1), xi_min(2,k+1), xi_max(3,k+1), xi_max(4,k+1);...
                        xi_min(1,k+1), xi_min(2,k+1), xi_max(3,k+1), xi_max(4,k+1);...
                        xi_min(1,k+1), xi_max(2,k+1), xi_max(3,k+1), xi_max(4,k+1);...
                        xi_max(1,k+1), xi_max(2,k+1), xi_max(3,k+1), xi_min(4,k+1);...
                        xi_max(1,k+1), xi_min(2,k+1), xi_max(3,k+1), xi_min(4,k+1);...
                        xi_min(1,k+1), xi_min(2,k+1), xi_max(3,k+1), xi_min(4,k+1);...
                        xi_min(1,k+1), xi_max(2,k+1), xi_max(3,k+1), xi_min(4,k+1);...
                        xi_max(1,k+1), xi_max(2,k+1), xi_min(3,k+1), xi_min(4,k+1);...
                        xi_max(1,k+1), xi_min(2,k+1), xi_min(3,k+1), xi_min(4,k+1);...
                        xi_min(1,k+1), xi_min(2,k+1), xi_min(3,k+1), xi_min(4,k+1);...
                        xi_min(1,k+1), xi_max(2,k+1), xi_min(3,k+1), xi_min(4,k+1);...
                        xi_max(1,k+1), xi_max(2,k+1), xi_min(3,k+1), xi_max(4,k+1);...
                        xi_max(1,k+1), xi_min(2,k+1), xi_min(3,k+1), xi_max(4,k+1);...
                        xi_min(1,k+1), xi_min(2,k+1), xi_min(3,k+1), xi_max(4,k+1);...
                        xi_min(1,k+1), xi_max(2,k+1), xi_min(3,k+1), xi_max(4,k+1);...
                        ];
                else
                    error('Sorry! Now we only accept the systems with 2 or 4 dimensions. Please make sure sys.nx = 2 or 4')
                end
                Uc_vertex = [ui_max(1,k); ui_min(1,k)];
                Xc_bar{k+1} = Polyhedron(Xc_vertex);
                Uc_bar{k} = Polyhedron(Uc_vertex);
            end
        end
        
        function Xf_bar = nominal_TerminalSet(obj,Xf)
            [F, ~, ~] = convert_Poly2Mat(Xf,obj.Uc);
            Hf = F(1:end-2,:); hf = ones(size(Hf,1),1);
            mu = obj.mu_e(:,end); Sigma = obj.Sigma_e{end};
            
            % X is the random variable with N(mu,Sigma)
            for j = 1:size(Hf,1) % number of state constraints
                mu_X(j,1) = Hf(j,:)*mu-hf(j,1);
                std_X(j,1) = sqrt(Hf(j,:)*Sigma*Hf(j,:)');
            end
            
            tol = 0.05; % probability of satisfaction is 1-epsilon
            z_x = norminv(1-tol);
            eta_f = -(mu_X + z_x*std_X);
            Xf_bar = Polyhedron('A', Hf, 'b', eta_f, 'Ae', [], 'be', []);
        end
        
    end
end