classdef ModelPredictiveControl < handle
    
    properties (SetAccess = private)
        sys % linear sys
        optcon; % optimal contol solver object
        Xc
        Uc
        Xc_bar
        Uc_bar
%         XMPI
        solution_cache
    end
    
    methods (Access = public)
        
        function obj = ModelPredictiveControl(sys, Xc, Uc, Xc_bar, Uc_bar, N)
            obj.sys = sys;
            obj.Xc_bar = Xc_bar; % time-varying state constraint set
            obj.Xc = Xc; % static and original state constraint set
            obj.Uc_bar = Uc_bar; % time-varying control constraint set
            obj.Uc = Uc; % static and original control constraint set
            obj.optcon = OptimalControler(sys, Xc_bar, Uc_bar, N);
%             obj.XMPI = obj.sys.compute_MPIset(Xc_bar{end}, Uc_bar{end});
            obj.solution_cache = [];
        end

        function u_next = solve(obj, x_init)
            obj.optcon.add_initial_eq_constraint(x_init);
            [x_nominal_seq, u_nominal_seq] = obj.optcon.solve();
            if isempty(u_nominal_seq) || isempty(x_nominal_seq)
                infeasible = true;
            else
                infeasible = false;
            end
            
            if infeasible
                u_next = []; x_nom = []; u_nominal = [];
            else
                obj.solution_cache = struct(...
                    'x_init', x_init', 'x_nominal_seq', x_nominal_seq, 'u_nominal_seq', u_nominal_seq);
                u_next = u_nominal_seq(:, 1);
            end
        end

        function show_prediction(obj,i,mu,Sigma,CI_Z,Xf_bar)
            assert(~isempty(obj.solution_cache), 'can be used only after solved');
            if iscell(obj.Xc_bar)
                Graphics.show_convex(obj.Xc,'Plot','yellow','FaceAlpha',0.1);
%                 Graphics.show_convex(obj.Xc_bar{i},'Plot',[0.8500 0.3250 0.0980],'FaceAlpha',0.4);
                
                Xmpi_robust = obj.sys.compute_MPIset(obj.Xc, obj.Uc);
                Graphics.show_convex(Xmpi_robust,'Plot',[0.75, 0.75, 0.75]); % gray
                Graphics.show_convex(Xf_bar,'Plot',[0.5, 0.5, 0.5]); % gray
                
%                 Xmpi_robust_bar = obj.sys.compute_MPIset(obj.Xc_bar{i}, obj.Uc_bar{i});
%                 Graphics.show_convex(Xmpi_robust_bar,'Plot',[0.5, 0.5, 0.5]);
                
                x_init = obj.solution_cache.x_init;
                

                x_nominal_seq = obj.solution_cache.x_nominal_seq;
                Graphics.show_trajectory(x_nominal_seq,'bs-','LineWidth',1.5,'MarkerSize',10);
                
                for j = 2:size(x_nominal_seq,2)
                    [C1,C2,C3] = find_multidim_contour(mu(:,j-1)',Sigma{j-1},CI_Z);
                    fill(x_nominal_seq(1,j)+C1(1,2:end),x_nominal_seq(2,j)+C1(2,2:end), 'g', 'FaceAlpha', 0.1, 'LineStyle','--', 'LineWidth', 0.5);
                    fill(x_nominal_seq(1,j)+C2(1,2:end),x_nominal_seq(2,j)+C2(2,2:end), 'g', 'FaceAlpha', 0.1, 'LineStyle','--', 'LineWidth', 0.5);
                    fill(x_nominal_seq(1,j)+C3(1,2:end),x_nominal_seq(2,j)+C3(2,2:end), 'g', 'FaceAlpha', 0.1, 'LineStyle','--', 'LineWidth', 0.5);
                end
                scatter(x_init(1), x_init(2), 50, 'ks', 'filled');
%                 leg = legend({'$X_c$', '$X_f$', 'nominal traj.', 'prob tubes', 'current state'}, 'Location', 'southeast');
            else
                Graphics.show_convex(obj.Xc_bar,'Plot','yellow','FaceAlpha',0.1);
            
                Xmpi_robust = obj.sys.compute_MPIset(obj.Xc_bar, obj.Uc_bar);
                Graphics.show_convex(Xmpi_robust,'Plot',[0.75, 0.75, 0.75]); % gray

                x_init = obj.solution_cache.x_init;
                scatter(x_init(1), x_init(2), 50, 'ks', 'filled');

                x_nominal_seq = obj.solution_cache.x_nominal_seq;
                Graphics.show_trajectory(x_nominal_seq, 'b.-','LineWidth',1.5,'MarkerSize',10);
                
%                 leg = legend({'$X_c$', '$X_f$', 'current state', 'nominal traj.'}, 'position', [0.5 0.15 0.1 0.2]);
            end

%             set(leg, 'Interpreter', 'latex');
        end
    end
    
end
