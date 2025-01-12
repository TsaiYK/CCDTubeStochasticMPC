function [prob_vio_state,prob_vio_ctrl] = PlotGenerator(System,Constraint,...
    mpc,x_nom,u_nom,x,u,x_for_plot,u_for_plot,i,n,savedir_name,PlotActualTraj)
if nargin<13
    PlotActualTraj = false;
    if nargin<12
        savedir_name = pwd;
        if nargin<11
            n = 1e2; % number of samples
        end
    end
end

%% Initialize trajectories
k_for_u = (i-1):(i+System.N_horizon-2);
k_for_x = (i-1):(i+System.N_horizon-1);

%% Initialize figures
figure(1)
for k = 1:System.nx
    subplot(System.nx,1,k); hold on
end
set(gcf,'Position',[680    42   560   954])
figure(2); hold on

%% Plot actual trajectories
num_plot = 100;
rand_plot_index = randi([1,n],1,num_plot);
if PlotActualTraj
    for k = 1:n
        if n<=num_plot % if sample size is not too big, plot all points
            figure(1)
            subplot(2,1,1)
            plot(k_for_x,x{k}(1,:),'--','LineWidth',0.5, 'Color',[0,0,1,0.2])
            subplot(2,1,2)
            plot(k_for_x,x{k}(2,:),'--','LineWidth',0.5, 'Color',[0,0,1,0.2])

            figure(2)
            plot(k_for_u,u{k},'--','LineWidth',0.5, 'Color','#FFAFAF')
        else % else plot 1e3 of all the points
            if sum(k==rand_plot_index)~=0
                figure(1)
                subplot(2,1,1)
                plot(k_for_x,x{k}(1,:),'--','LineWidth',0.5, 'Color',[0,0,1,0.2])
                subplot(2,1,2)
                plot(k_for_x,x{k}(2,:),'--','LineWidth',0.5, 'Color',[0,0,1,0.2])

                figure(2)
                plot(k_for_u,u{k},'--','LineWidth',0.5, 'Color','#FFAFAF')
            end
        end
    end
end

%% Plot the nominal trajectories
figure(1)
for j = 1:System.nx
    subplot(System.nx,1,j); hold on
    % state limits (lower and upper bounds)
    plot([(i-1), System.N_horizon+10],[System.xmin(j),System.xmin(j)],'k--','LineWidth',2)
    plot([(i-1), System.N_horizon+10],[System.xmax(j),System.xmax(j)],'k--','LineWidth',2)
    % fill the feasible region of state
    fill([0:System.N fliplr(0:System.N)],...
        [ones(1,System.N+1)*System.xmin(j),ones(1,System.N+1)*System.xmax(j)],...
        'y', 'FaceAlpha', 0.1,'EdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1)
    plot(k_for_x,x_nom(j,:),'b--','LineWidth',1,'MarkerSize',5);
    xlim([0,(System.N)])
    ylim([System.xmin(j),System.xmax(j)]*1.1)
    ylabel(strcat('$\xi^{(',num2str(j),')}_k$'))
end
subplot(System.nx,1,System.nx); xlabel('$k$')

figure(2)
plot([(i-1), System.N_horizon+10],[System.umin, System.umin],'k--','LineWidth',2)
plot([(i-1), System.N_horizon+10],[System.umax, System.umax],'k--','LineWidth',2)
fill([0:System.N fliplr(0:System.N)],...
        [ones(1,System.N+1)*System.umin,ones(1,System.N+1)*System.umax],...
    [0.8500 0.3250 0.0980], 'FaceAlpha', 0.1,...
    'EdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1)
p2 = plot(k_for_u,u_nom,'k--','LineWidth',1,'MarkerSize',5,'Color','#CC0000');
ylim([System.umin,System.umax])
xlim([0,(System.N)])
xlabel('$k$')
ylabel('$u_k$')

% Plot the 90%, 95%, 99% CIs
CI_Z = [1.6449, 1.96, 2.5758]; % Z values for the 90%, 95%, 99% CIs
for j = 2:(System.N_horizon+1)
    [C{j,1},C{j,2},C{j,3}] = find_multidim_contour(Constraint.mu_e(:,j-1)',...
        Constraint.Sigma_e{j-1},CI_Z);
end
for num = 1:3
    min_state_CI{num} = zeros(System.nx,System.N_horizon+1);
    min_state_CI{num}(:,1) = x_nom(:,1);
    max_state_CI{num} = zeros(System.nx,System.N_horizon+1);
    max_state_CI{num}(:,1) = x_nom(:,1);
    for j = 2:(System.N_horizon+1)
        min_state_CI{num}(:,j) = x_nom(:,j) + min(C{j,num},[],2);
        max_state_CI{num}(:,j) = x_nom(:,j) + max(C{j,num},[],2);
    end
    
    min_control_CI{num} = zeros(System.nu,System.N_horizon);
    min_control_CI{num}(:,1) = u_nom(:,1);
    max_control_CI{num} = zeros(System.nu,System.N_horizon);
    max_control_CI{num}(:,1) = u_nom(:,1);
    for j = 1:System.N_horizon
        min_control_CI{num}(:,j) = u_nom(:,j) - CI_Z(num)*sqrt(Constraint.Sigma_u(j));
        max_control_CI{num}(:,j) = u_nom(:,j) + CI_Z(num)*sqrt(Constraint.Sigma_u(j));
    end
end

figure(1)
for j = 1:System.nx
    subplot(System.nx,1,j)
    for num = 1:3
        fill([k_for_x fliplr(k_for_x)],...
            [min_state_CI{num}(j,:) fliplr(max_state_CI{num}(j,:))], 'b', 'FaceAlpha', 0.2,...
            'EdgeColor', 'none', 'LineStyle', '--', 'LineWidth', 0.5)
    end
end

figure(2)
for num = 1:3
    fill([k_for_u fliplr(k_for_u)],...
        [min_control_CI{num} fliplr(max_control_CI{num})], 'r', 'FaceAlpha', 0.2 ,...
        'EdgeColor', 'none', 'LineStyle', '--', 'LineWidth', 0.5)
end

vio_count_state = zeros(size(System.H,1),System.N_horizon+1);
vio_count_control = zeros(size(System.G,1),System.N_horizon);
tic
for k = 1:n
    if PlotActualTraj
        if n<=num_plot % if sample size is not too big, plot all points
            scatter(x{k}(1, :), x{k}(2, :), 10, 'MarkerFaceAlpha',0.2,...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor','b','HandleVisibility','off')
        else % else plot 1e3 of all the points
            if sum(k==rand_plot_index)~=0
                scatter(x{k}(1, :), x{k}(2, :), 10, 'MarkerFaceAlpha',0.2,...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor','b','HandleVisibility','off')
            end
        end
    end
    g1 = System.H*x{k}-System.h;
    boolean_g1 = g1>0;
    vio_count_state(boolean_g1) = vio_count_state(boolean_g1)+1;
    g2 = System.G*u{k}-System.g;
    boolean_g2 = g2>0;
    vio_count_control(boolean_g2) = vio_count_control(boolean_g2)+1;
end
toc

prob_vio_state = vio_count_state/n;
prob_vio_ctrl = vio_count_control/n;

%% Past trajectories
figure(1)
for j = 1:System.nx
    subplot(System.nx,1,j)
    plot(0:(i-1),x_for_plot(j,:),'b.-','MarkerSize',10)
    scatter(i-1,x_for_plot(j,end),50,'ks','filled');
end
filename = strcat(savedir_name, 'tmpc_state_seq', number2string(i), '.png');
saveas(gcf, char(filename)); % removing this line makes the code much faster
close(gcf);

figure(2)
if i~=1
    plot((0:(i-1)),u_for_plot,'r.-','MarkerSize',10)
end
p1 = scatter(i-1,u_for_plot(end),50,'ks','filled');
filename = strcat(savedir_name, 'tmpc_ctrl_seq', number2string(i), '.png');
legend([p1,p2],{'$u_k$','$\bar{\mathbf{U}}$'})
saveas(gcf, char(filename)); % removing this line makes the code much faster
close(gcf);

if System.nx==2
    % Plot trajectory on a phase plane
    figure
    mpc.show_prediction(i,Constraint.mu_e,Constraint.Sigma_e,CI_Z,Constraint.Xf_bar);
    axis([-11,6,-6,3])
    xlabel('$\xi^{(1)}$')
    ylabel('$\xi^{(2)}$')
    filename = strcat(savedir_name, 'tmpc_plane_seq', number2string(i), '.png');
    saveas(gcf, char(filename)); % removing this line makes the code much faster
    close(gcf);
end

end