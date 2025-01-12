function [x_opt,fval,M] = runobjconstr(xp,SaveLoc_dir,DesignStrategy,i,W_constr_handling)

xLast = []; % Last place computeall was called
myf1 = []; % Use for objective at xLast
myf2 = []; % Use for objective at xLast
myf3 = []; % Use for constraint at xLast

A = []; b = []; Aeq = []; beq = [];
TolCon = 1e-6; classif_err_allowance = 0.1;

%% P3GA Setting
fun_all = @(x) [myfun1(x) myfun2(x)];    % objective function
par = [];                          % the parameter function is in the first column of fun1and2
dom = [1,2];                          % the dominator function is in the second column of fun1and2
if strcmp(DesignStrategy,'CCD')
    K_lb = [-1.2583,-2.2467];
    K_ub = [-0.0948,-0.2351];
    lb = [-1,-1,K_lb];
    ub = [1,1,K_ub];
elseif strcmp(DesignStrategy,'plantdesign')
    K_lqr = [-0.6167,-1.2703];
    % For Plant Design:
    lb = [-1,-1];
    ub = [1,1];
elseif strcmp(DesignStrategy,'controldesign') || strcmp(DesignStrategy,'seqdesign')
    K_lb = [-1.2583,-2.2467];
    K_ub = [-0.0948,-0.2351];
    lb = [K_lb];
    ub = [K_ub];
else
    error('Unrecognizable design strategy!')
end

nvars = length(lb);                  % the number of design variables
Generations = 3;                    % maximum number of generations
PopulationSize = 5;                % maximum populations
nonlcon = @(x) myfun3(x);
% The constraint is set empty here because we assign NaN to both objectives
% if the design is infeasible (including no infeasible solution for nominal
% MPC, no disturbance set Z exists, and AK(=A+BK) is not exponentially stable)

if ~W_constr_handling
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true);
    options.SaveLoc = SaveLoc_dir;
else
    % Initialize population
    [X_new, feasible, initial, initialfeasible, cEval] = Initialize_LHS_GP(lb,ub,...
            A,b,Aeq,beq,nonlcon,PopulationSize,nvars,'DontPlot',classif_err_allowance,TolCon,[]);
    % Define options
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true,'InitialPopulation',X_new,...
        'CrossoverFraction',0.8,'MutationFraction',0.05);
    options.initial = initial;
    options.initialfeasible = initialfeasible;
    options.GenerationPlots = 1;
    options.classif_err_allowance = classif_err_allowance;
    options.cEval = cEval;
    options.TolCon = TolCon;
    options.SaveLoc = SaveLoc_dir;
end
options.Log = true;
options.GenerationData = true;
options.GenerationDataIncrement = 2;
options.Dominance = 'hch'; % determine if you want to use HCH-based dominance
options.hvi_options.phi = 90; % Hypercone angle
options.ref_bnds = [200,600;0,0.5];

%% Run P3GA
[x_opt,fval,M] = p3ga(fun_all,dom,par,nvars,[],[],[],[],lb,ub,nonlcon,options);
figure(gcf); 
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_phvi.fig'))
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_phvi.png'))
close(gcf)

figure(gcf); 
xlabel('$V_{N_s}$','Interpreter','Latex');
ylabel('$||$eig($\Sigma_{\mathbf{e}})||_2$','Interpreter','Latex');
zlabel('$u_{max}$','Interpreter','Latex')
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_results.fig'))
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_results.png'))
close(gcf)

if strcmp(DesignStrategy,'seqdesign')
    save(strcat(options.SaveLoc,'\p3ga_results',num2str(i),'.mat'))
else
    save(strcat(options.SaveLoc,'\p3ga_results.mat'))
end

%% Private functions for objective
    function [f1] = myfun1(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            if strcmp(DesignStrategy,'plantdesign')
                % For plant design
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x;
                K = K_lqr;
            elseif strcmp(DesignStrategy,'controldesign') || strcmp(DesignStrategy,'seqdesign')
                % For Control Design
                K = x(1:2);
            else
                % For CCD
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x(1);
                K = x([3,4]);
            end
            [J,Size_Sigma_e] = func_tubeSMPC(xp,K,'2D_numerical');
            myf1 = J;
            myf2 = Size_Sigma_e;
            myf3 = isnan(J);
            xLast = x;
        end
        % Now compute objective function
        f1 = myf1;
    end


    function [f2] = myfun2(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            if strcmp(DesignStrategy,'plantdesign')
                % For plant design
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x;
                K = K_lqr;
            elseif strcmp(DesignStrategy,'controldesign') || strcmp(DesignStrategy,'seqdesign')
                % For Control Design
                K = x(1:2);
            else
                % For CCD
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x(1);
                K = x([3,4]);
            end
            [J,Size_Sigma_e] = func_tubeSMPC(xp,K,'2D_numerical');
            myf1 = J;
            myf2 = Size_Sigma_e;
            myf3 = isnan(J);
            xLast = x;
        end
        % Now compute objective function
        f2 = myf2;
    end

    function [c,ceq] = myfun3(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            if strcmp(DesignStrategy,'plantdesign')
                % For plant design
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x;
                K = K_lqr;
            elseif strcmp(DesignStrategy,'controldesign') || strcmp(DesignStrategy,'seqdesign')
                % For Control Design
                K = x(1:2);
            else
                % For CCD
                xp = 2.25*x(1)*sin(x(2)*pi)+2.75;
%                 xp = x(1);
                K = x([3,4]);
            end
            [J,Size_Sigma_e] = func_tubeSMPC(xp,K,'2D_numerical');
            myf1 = J;
            myf2 = Size_Sigma_e;
            myf3 = isnan(J);
            xLast = x;
        end
        % Now compute constraint function
        if myf3 % infeasible design
            c = 1;
        else
            c = -1;
        end
        ceq = [];
    end
    
end
