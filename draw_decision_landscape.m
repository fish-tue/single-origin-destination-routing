%% DRAW_DECISION_LANDSCAPE - draw decision landscape
% Description:
%   Draw solution landcape according to Theorems III.1 and III.2 in [1].
% Outputs: none
% Assumptions and limitations: 
%   - k_ref = 0
% Other m-files required:
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

%% Draw solution to individual users' problem
% Assumption k_ref = 0
k_ref = 0;
% Range of k
k = (-p(1):.1:k_ij(1,1,k_ref,p,T)+p(1))';
% Compute gamma_j, j = 1,...,n with gamma_n_arcs_individual
gamma = gamma_n_arcs_individual(d_star,T,p,k,k_ref,s_min,s_bar,s_max);
% Decision Evolution
figure('Position',4*[0 0 192 144]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
k_inf = max([0,k_ij(n,n,k_ref,p,T)]);
feasible_mask = k >= k_inf;
unfeasible_mask = k <= k_inf;
for j = 1:n
    if j == 1
        % Decision 1
        aux_x = k(feasible_mask);
        aux_y = (s_max/s_bar)*ones(size(k(feasible_mask),1),1);
        area(aux_x,aux_y,'LineWidth',2);
    else
        % Decision 2-n
        aux_x = k(feasible_mask);
        aux_y = gamma(j-1,feasible_mask);
        area(aux_x,aux_y,'LineWidth',2);
    end
end
% Unfeasible
aux_x = k(unfeasible_mask);
aux_y = (s_max/s_bar)*ones(size(k(unfeasible_mask),1),1);
area(aux_x,aux_y,'LineWidth',2,'FaceColor',[0 0 0],'FaceAlpha',0.3);
% Lines for a pretty plot
plot(k_inf*ones(2,1),[(s_min/s_bar);(s_max/s_bar)],'LineWidth',2,...
    'Color','black')
legend({' Arc 1',' Arc 2',' Arc 3',' Arc 4',' Arc 5',' Unfeasible',},...
    'Interpreter','latex','Location','northeast');
ylabel('$s/\bar{s}$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
xlim([k(1) k(end)]);
ylim([(s_min/s_bar) (s_max/s_bar)]);
hold off;
% Save figure to .fig and .eps formats
savefig('./fig/decision_landscape.fig')
set(gcf,'renderer','Painters')
saveas(gca,'./fig/decision_landscape.eps','epsc');
clear;

%% Auxiliary functions
% Karma thresholds for unitary decisions, i.e., [y_bar]_j = 1 for some j
function k_ij = k_ij(i,j,k_ref,p,T)
    k_ij = k_ref+p(i)+T*p(j);
end
