%% simulation - simulate incentive mechanism
% Description:
%   Simulation of the daily Nash equilibrium and the decisions of each user
% Outputs:
%   simulation.mat (simulation output)
% Assumptions and limitations:
%   - Sensitivity distribution is uniform
%   - Initial karma distribution is dicrete uniform with support 
%   {25p_1,25p_1,+1,...,50p_1}
%   - k_ref is a discrete uniform distribution with support {k_ref \in N : 
%     k_ref = 0 \lor k_ref = p_j, j = 1,...,n}
% Other m-files required:
%   - n_arcs_individual.m
%   - gamma_n_arcs_individual.m
% MAT-files required:
%   - network.mat (generated with generate_network.m)
%   - prices.mat (computed with n_arcs_pricing.mat)
% Toolboxes required: none
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W. P. M. H. Heemels, and M. Salazar, “Urgency-aware 
%   optimal routing in repeated games through artificial currencies” 
%   [not published yet]

%% Initialization
clear;
% Load illustrative network
load('network.mat','n','M','s_min','s_bar','s_max','P_home','P_go','T',...
    'alpha','beta','d0','kappa','d','c0','c','x_star','d_star');
% Load prices
load ('prices.mat','p');

%% Parameters
% Nash iterations
it_nash_max = 50;
it_nash_epsl = 1e-3;
% Simulation
T_sim = 200; 
k_init_max = 50*max(p);

%% Simulation
% Simulation - Initialization
k = zeros(M,T_sim+1); % karma 
s = zeros(M,T_sim+1); % sensitivity
y = zeros(M,T_sim+1);
x = zeros(n,T_sim+1); % aggregate choices
it_nash = zeros(T_sim+1,1); % number of iterations for the nash equilibrium
% Karma initialization
k(:,1) = round(unifrnd(k_init_max/2,k_init_max,M,1)); % uniform initial karma distribution
% k_ref initialization
p_p = round([0; p(p>0)]);
k_ref_idx = unidrnd(length(p_p),M,1);
k_ref = p_p(k_ref_idx);
% Simulation loop
for t = 1:T_sim+1
    % Pick daily sensitivities
    s(:,t) = unifrnd(0,s_max,M,1);
    % Agents decisions
    y_go = binornd(1,P_go,M,1); % mask of travelling agents    
    % Init Nash equilibrium iterations
    if t > 1      
        x(:,t) = x(:,t-1);
    else
        x(:,t) = x_star;
    end
    % Until convergence of Nash equilibrium is reached
    I = eye(n);
    y(:,t) = ones(M,1);
    while true
        % Previous iteration's flows
        x_prev = x(:,t);
        % Agents decisions
        for i = 1:M
            if ~y_go(i)
                y(i,t) = 0;
                continue; 
            end % non-traveling agent
            y(i,t) = n_arcs_individual(...
               d(x(:,t)-I(:,y(i,t))/M+ones(n,1)/M),T,p,k(i,t),k_ref(i),s(i,t),s_min,s_bar,s_max);
            % Aggregate behaviour
            for j = 1:n
                x(j,t) = sum(y(:,t)==j)/M;
            end
        end
        % Catch infeasibility
        if sum(isnan(y(:,t)))
            error("Caught infeasibility!");
        end
        % Catch Nash equilibrium iterations not converging
        it_nash(t) = it_nash(t)+1;
        if norm(x_prev-x(:,t))<1e-3/M, break; end
        if it_nash(t) > it_nash_max
            warning("Nash iterations did not converge.");
            break; 
        end
    end
    % Karma dynamics
    if t < T_sim+1
        for i = 1:M
            if ~y_go(i)
                k(i,t+1) = k(i,t);
            else
                k(i,t+1) = k(i,t)-p(y(i,t));
            end
        end
    end  
end

%% Plot simulation evolution
% Decision Evolution
figure('Position',4*[0 0 192 144]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
for j = 1:n
    aux_x = (0:T_sim)';
    aux_y = sum(x(j:n,:),1)';
    [aux_x,aux_y] = stairs_vector(aux_x,aux_y);
    area(aux_x,aux_y,'LineWidth',2);
end
for j = 1:n
    plot([0 T_sim],[sum(x_star(j:n)) sum(x_star(j:n))],'--','Color','black','LineWidth',3);
end
legend({' Arc 1',' Arc 2',' Arc 3',' Arc 4',' Arc 5',' $\,x^\star$'},'Location','southeast','Interpreter','latex');
ylabel('$\mathbf{x}^\mathrm{NE}(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([0 1]);
hold off;
% Save figure to .fig and .eps formats
savefig('./fig/decision.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/decision.eps','epsc');

% Karma levels
k_mean = mean(k,1);
k_max = max(k);
k_min = min(k);
k_var = sqrt(var(k));
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
aux_x2 = [0:T_sim fliplr(0:T_sim)]';
aux_y2 = [k_max fliplr(k_min)]';
[aux_x2,aux_y2] = stairs_vector(aux_x2,aux_y2);
fill(aux_x2,aux_y2*1e-2,'k',...
    'LineWidth',2,'FaceColor','black','FaceAlpha',0.2,'EdgeAlpha',0);
aux_x1 = [0:T_sim fliplr(0:T_sim)]';
aux_y1 = [k_mean+k_var fliplr(k_mean-k_var)]';
aux_y1(aux_y1<0) = 0;
[aux_x1,aux_y1] = stairs_vector(aux_x1,aux_y1);
fill(aux_x1,aux_y1*1e-2,'k',...
    'LineWidth',2,'FaceColor','black','FaceAlpha',0.4,'EdgeAlpha',0);
stairs(0:T_sim,k_mean*1e-2,'LineWidth',2,'Color','black');
stairs(aux_x1,aux_y1*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
stairs(aux_x2,aux_y2*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
legend({' $\max_i$/$\min_i$ $\{k^i(t)\}$',' $\hat{k}(t)\pm\sigma_k(t)$',' $\hat{k}(t)$'},...
    'Location','northeast','Interpreter','latex');
%ylabel('$k^i(t), i = 1,\ldots,M$','Interpreter','latex');
ylabel('Karma level $\times 10^{-2}$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
hold off;
% Save figure to .fig and .eps formats
savefig('./fig/karma.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/karma.eps','epsc');

% System's cost
cost_soc = zeros(T_sim+1,1);
cost_soc_rel_opt = zeros(T_sim+1,1);
cost_soc_opt = c(x_star);
for t = 1:T_sim+1
    cost_soc(t) = c(x(:,t));
    cost_soc_rel_opt(t) = (cost_soc(t)-cost_soc_opt)/cost_soc_opt;
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
stairs(0:T_sim,cost_soc_rel_opt*100,'LineWidth',2.5,'Color','black');
ylabel('$\Delta$ Societal cost $(\%)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([-1 30]);
hold off;
% Save figure to .fig and .eps formats
savefig('./fig/cost.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/cost.eps','epsc');

% Sensitivity to discomfort
delta_s = 100*mean(s-s_bar)/s_bar;
delta_d = zeros(T_sim+1,1);
I = eye(n);
for t = 1:T_sim+1
    aux1 = 0;
    aux2 = 0;
    for i = 1:M
        if ~y(i,t), continue; end
        aux1 = aux1 + s(i,t)*d(x(:,t))'*I(:,y(i,t)) - s_bar.*d(x(:,t))'*I(:,y(i,t));
        aux2 = aux2 + s_bar.*d(x(:,t))'*I(:,y(i,t));
    end
    delta_d(t) = 100*aux1/aux2;
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex');
stairs(0:T_sim,delta_s,':','LineWidth',2.5,'Color','black');
stairs(0:T_sim,delta_d,'LineWidth',2.5,'Color','black');
legend({'$\Delta \bar{s}(t)$','$\Delta \bar{d}(t)$'},'Interpreter','latex','Location','southwest');
ylabel('Sensitivity, discomfort [$\%$]','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
hold off;
% Save figure to .fig and .eps formats
savefig('./fig/sensitivity.fig');
set(gcf,'renderer','Painters')
saveas(gca,'./fig/sensitivity.eps','epsc');

%% Save simulation results
save('simulation.mat')
clear;

%% Auxiliary function
% To user the equivalent of stairs with area and fill functions
function [x,y] = stairs_vector(x,y)
    x = [x';x'];
    y = [y';y'];
    y = y(:);
    x = x([2:end end])';
end
