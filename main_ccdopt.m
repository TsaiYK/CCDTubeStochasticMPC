clear
clc
close all

%% Add paths
addpath('../src/')
addpath('../src/utils/')
W_constr_handling = true;
if W_constr_handling
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
else
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
end
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc'))
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc\resources\cprnd')) % for some reason this doesn't work sometimes
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\EPO'));
addpath(genpath('C:\Users\yktsai0121\OneDrive - Texas A&M University\2022Fall\IEEE_Evo_Comp_2023\matlab_codes'))

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1.5);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',14)

%% Select a Design Strategy
DesignStrategy = 'seqdesign';
SystemType = '2D_numerical';

mkdir('results')

if strcmp(DesignStrategy,'CCD')
    SaveLoc_dir = 'results\ccd_results';
    xp = [];
elseif strcmp(DesignStrategy,'plantdesign')
    SaveLoc_dir = 'results\plantdesign_results';
    xp = [];
elseif strcmp(DesignStrategy,'controldesign')
    SaveLoc_dir = 'results\controldesign_results';
    switch SystemType
        case '2D_numerical'
            xp = 1; % if control design only, the plant design variable is fixed
        case '4D_satellite'
            xp = cosd(39.14); % if control design only, the plant design variable is fixed
    end
    
elseif strcmp(DesignStrategy,'seqdesign')
    % if plant design is done and we are going to do control design based on
    % all the designed plants
    SaveLoc_dir = 'results\seqdesign_results';
    PlantOptDesign = load('results\plantdesign_results\xp_opt.mat'); % from the plant design optimization results
    xp_opt = PlantOptDesign.xp_opt;
    if strcmp(SystemType,'2D_numerical')
%         theta_p_opt = 2.25*xp_opt(:,1).*sin(xp_opt(:,2)*pi)+2.75;
%         xp_opt = PlantOptDesign.xp_opt;
        rand_index = randi([1,size(xp_opt,1)],1,50);
        theta_p_opt = 2.25*xp_opt(rand_index,1).*sin(xp_opt(rand_index,2)*pi)+2.75;
    else
        theta_p_opt = xp_opt;
    end
end

if strcmp(DesignStrategy,'seqdesign')
    x_opt = []; fval = [];
    for i = 1:length(theta_p_opt)
        [K_opt{i},fval_iter{i},M{i}] = runobjconstr(theta_p_opt(i),SaveLoc_dir,DesignStrategy,i,W_constr_handling);
        x_opt_tmp = [repmat(theta_p_opt(i),size(K_opt{i},1),1),K_opt{i}];
        x_opt = [x_opt;x_opt_tmp];
        fval = [fval;fval_iter{i}]; 
    end
else
    [x_opt,fval,M] = runobjconstr_p3ga(xp,SaveLoc_dir,DesignStrategy,[],W_constr_handling);
end

%% Approximate PPD
f1 = fval(:,1);
f2 = fval(:,2);

B = [f1,f2];
x_fea = x_opt;
x_fea(isnan(B(:,1)),:) = [];
B(isnan(B(:,1)),:) = [];
[A, varargout] = prtp(B);
f1_ndom = A(:,1);
f2_ndom = A(:,2);
x_ndom = x_fea(varargout,:);

%% Plot PF
figure(1)
hold on
% p2 = plot(f1_PPD,f2_PPD,'ko');
xlabel('$V_{N_s}$','Interpreter','Latex');
ylabel('$||$eig($\Sigma_{\mathbf{e}})||$','Interpreter','Latex');

[f1_PPD_sort, index_sort] = sort(f1_ndom);
f2_PPD_sort = f2_ndom(index_sort);
x_PPD_sort = x_ndom(index_sort,:);

if strcmp(DesignStrategy,'CCD')
    plot(f1_PPD_sort,f2_PPD_sort,'b-*'); % CCD
    title('Design Strategy: CCD')
elseif strcmp(DesignStrategy,'plantdesign')
    plot(f1_PPD_sort,f2_PPD_sort,'-^','Color',[0.3010 0.7450 0.9330]); % only plant design
    title('Design Strategy: Plant Design')
elseif strcmp(DesignStrategy,'controldesign')
    plot(f1_PPD_sort,f2_PPD_sort,'k-o'); % control design by fixing plant
    title('Design Strategy: Control Design')
elseif strcmp(DesignStrategy,'seqdesign')
    plot(f1_PPD_sort,f2_PPD_sort,'-^','Color',[0.4660 0.6740 0.1880]); % seq design
    title('Design Strategy: Sequential Design')
end

%% Save data
save(strcat(SaveLoc_dir,'\p3ga_results_final.mat'))
if strcmp(DesignStrategy,'plantdesign')
    xp_opt = x_opt;
    save(strcat(SaveLoc_dir,'\xp_opt.mat'),'xp_opt')
end
