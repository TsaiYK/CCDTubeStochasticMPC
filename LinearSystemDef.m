classdef LinearSystemDef < handle
    properties (SetAccess = public)
        A         % matrix A for system dynamics (inertia)
        B         % matrix B for system dynamics (controllability)
        Gd        % matrix Gd for system dynamics (disturbance)
        Q         % weighting matrix for state
        R         % weighting matrix for input
        mysys     % linear system defined by A, B, Q, R, K
        K_lqr     % LQR gain
        A_K       % closed-loop system matrix A_K:=A+B*K
        iniCon    % initial condition/state
        N         % simulation period
        N_horizon % horizon
        nx
        nu
        xmin
        xmax
        umin
        umax
        H
        h
        G
        g
    end
    
    %% Public functions
    methods (Access = public)
        % 2-D Numerical example
        function obj = numerical_example(obj, xp, K, control_limit)
            %% System parameters
            obj.A = [1 1; 0 1];
            obj.B = [0.5; xp]; 
            obj.Q = diag([1, 1]);
            obj.R = 0.1;
            obj.K_lqr = -dlqr(obj.A,obj.B,obj.Q,obj.R); % please keep in mind that a "minus" should be included
            obj.iniCon = [-7; -2];
            obj.N = 20;
            obj.N_horizon = 20;

            % Define the linear system
            obj.mysys = LinearSystem(obj.A, obj.B, obj.Q, obj.R, K);

            % Define A_K = A+B*K
            obj.A_K = obj.A+obj.B*K;

            %% State and control constraints
            obj.xmin = [-10,-5];
            obj.xmax = [5,2];
            obj.umin = -control_limit;
            obj.umax = control_limit;
            [obj.H, obj.h, obj.G, obj.g, obj.nx, obj.nu] = obj.IneqConstrMatrix();
        end
        
        % 4-D satellite (roll-yaw orientation system)
        function obj = satellite_dynamics(obj, xp, K, control_limit)
            %% System parameters
            T = 0.1;
            Iy = 1; Ix = 3; Iz = 3.5;
            kx = (Ix-Iy)/Ix;
            ky = (Iz-Ix)/Iy;
            cos_alpha = xp(1);
            cos_beta = sqrt(1-cos_alpha^2);

            Ac = [0,1,0,0;...
                -4*kx,0,0,(kx-1);...
                0,0,0,1;...
                0,(ky+1),ky,0];
            Bc = [0;Iz/Ix*cos_alpha;0;Iz/Iy*cos_beta];

            funB = @(t) expm(Ac*t)*Bc;
            obj.A = expm(Ac*T);
            obj.B = integral(funB,0,T,'ArrayValued',true);
            obj.Q = zeros(4,4); obj.Q(1,1) = 1; obj.Q(3,3) = 1;
            obj.R = 1;
            obj.K_lqr = -dlqr(obj.A,obj.B,obj.Q,obj.R); % please keep in mind that a "minus" should be included
            obj.iniCon = [deg2rad(45);0.1;deg2rad(45);0.1];
            obj.N = 200;
            obj.N_horizon = 20;

            % Define the linear system
            obj.mysys = LinearSystem(obj.A, obj.B, obj.Q, obj.R, K);

            % Define A_K = A+B*K
            obj.A_K = obj.A+obj.B*K;

            %% State and control constraints
            obj.xmin = [-pi/2,-1.5,-pi/2,-0.75];
            obj.xmax = [pi/2,1.5,pi/2,0.75];
            obj.umin = -control_limit;
            obj.umax = control_limit;
            [obj.H, obj.h, obj.G, obj.g, obj.nx, obj.nu] = obj.IneqConstrMatrix();
        end
    end
    
    %% Private functions
    methods (Access = private)
        function [H,h,G,g,nx,nu] = IneqConstrMatrix(obj)
            % Transforming the lower and upper bounds to matrices
            nx = length(obj.xmin);
            nu = length(obj.umin);
            H = zeros(2*nx,nx); h = zeros(2*nx,1);
            for i = 1:nx
                H((i-1)*2+1,i) = 1;
                H(i*2,i) = -1;
                h((i-1)*2+1,1) = obj.xmax(i);
                h(i*2,1) = -obj.xmin(i);
            end
            G = zeros(2*nu,nu); g = zeros(2*nu,1);
            for i = 1:nu
                G((i-1)*2+1,i) = 1;
                G(i*2,i) = -1;
                g((i-1)*2+1,1) = obj.umax(i);
                g(i*2,1) = -obj.umin(i);
            end
        end
    end
end