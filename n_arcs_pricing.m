%% N_ARCS_PRICING - compute near-optimal arcs' prices
% Description:
%   Computes the near-optimal arc prices according to Section V.A in [1].
% Outputs:
%   prices.mat (n dimensional vector of prices, 'p')
% Assumptions and limitations:
%   - Sensitivity distribution is uniform
%   - k_ref is a discrete uniform distribution with support {k_ref \in N : 
%   k_ref = 0 \lor k_ref = p_j, j = 1,...,n}
%   - k_0 = p_1
% Other m-files required:
%   - gamma_n_arcs_individual.m
% MAT-files required:
%   - network.mat (generated with generate_network.m)
% Toolboxes required:
%   - Global Optimization Toolbox
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   05/02/2024 - Leonardo Pedroso
%       * Added final publication reference to [1] 
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W.P.M.H. Heemels and M. Salazar, "Urgency-Aware Routing
%   in Single Origin-Destination Itineraries Through Artificial 
%   Currencies," 2023 62nd IEEE Conference on Decision and Control (CDC), 
%   Singapore, Singapore, 2023, pp. 4142-4149, 
%   doi: 10.1109/CDC49753.2023.10383739.

%% Initialization
clear;
% Load illustrative network
load('network.mat','n','M','s_min','s_bar','s_max','P_home','P_go','T',...
    'alpha','beta','d0','kappa','d','c0','c','x_star','d_star');

%% Solve pricing problem - Genetic Alg. 
% Design parameters
% Bounds on price magnitude
ub_p = 100; 
lb_p = -ub_p;
x_star_quant_exp = 3;
% Quantization
x_star_quant = round((10^x_star_quant_exp)*x_star)/(10^x_star_quant_exp);
% Cost function
F = @(p) c(n_arcs_stationary_flows(d_star,T,p',s_max,s_bar,s_min,P_go,P_home));
% Linear inequality constraints
% p_1 > 0
% p_n < 0
% p_j > p_j+1 (p_j >= p_j+1 + 1)
c_ = zeros(n-1,1);
c_(1) = 1;
r = zeros(n,1);
r(1) = 1;
r(2) = -1;
s = zeros(n,1);
s(1) = -1;
s_ = zeros(n,1);
s_(end) = 1;
A = [s'; s_'; -toeplitz(c_,r)];
b = [-1;-1;-ones(n-1,1)];
% Linear equality constraints
% p^Tx_star = 0
Aeq = x_star_quant';
beq = 0;
% Price bounds
lb = lb_p;
ub = ub_p;
% Nonlinear inequality constraints
nolcon = [];
% Faster population generation: x_star perpendicular basis
x_star_pp_basis = zeros(n,n-1);
I = eye(n);
for i = 1:n-1
    for j = 1:i
        if j ~= 1
            aux = aux - x_star_pp_basis(:,j-1)*(x_star_pp_basis(:,j-1)'*I(:,i));
        else
            aux = I(:,i) - x_star_quant*(x_star_quant'*I(:,i)/(x_star_quant'*x_star_quant));
        end

    end
    x_star_pp_basis(:,i) = aux/sqrt(aux'*aux); 
end

% Hyperparameters
ga_population_size = 20;
ga_threshold = 1e-2;
ga_stall_it = 3;
% Stopping criterion
% - Cost below c(x_star)*(1+ ga_threshold)
% - Or constnat cost for ga_stall_it generations

% GA
fprintf("---------------------------------------------------------------------");
tic;
options = optimoptions('ga','Display','diagnose',...
    'CreationFcn', {@CustomCreationFcn, x_star_pp_basis, ub_p},...
    'PopulationSize',ga_population_size,'FitnessLimit',c(x_star)*(1+ga_threshold));
options.MaxStallGenerations = ga_stall_it;
p_opt_ga = ga(F,n,A,b,Aeq,beq,lb,ub,nolcon,1:n,options)';
p_opt_ga = round(p_opt_ga);
toc;
fprintf("---------------------------------------------------------------------\n");

% Output GA solution
fprintf("Genetic alg. prices:\t\t");
fprintf("%d ",p_opt_ga);
fprintf("\n");

%% Save pricing results
p = p_opt_ga;
save('prices.mat','p');
clear;

%% Auxiliary functions

% Aggregate stategy over distribution of k_ref
% (k_ref is a discrete uniform distribution with support {k_ref \in N : 
% k_ref = 0 \lor k_ref = p_j, j = 1,...,n})
function x_inf = n_arcs_stationary_flows(d,T,p,s_max,s_bar,s_min,P_go,P_home)
    % Support of theta_p distribution (distribution of k_ref values)
    k_ref = round([0;p(p>0)]);
    x_inf = zeros(length(p),length(k_ref));
    % Compute x_inf for each value in the support of theta_p
    for j = 1:length(k_ref)
        x_inf(:,j) = n_arcs_stationary_flows_kref(d,T,p,k_ref(j),s_max,s_bar,s_min,P_go,P_home);
    end
    % Expected aggregate (because theta_p is uniform)
    x_inf = mean(x_inf,2);
end

% Aggregate stategy for a particular value of k_ref
% (k_ref is a discrete uniform distribution with support {k_ref \in N : 
% k_ref = 0 \lor k_ref = p_j, j = 1,...,n})
% Assumptions on inputs:
%   - d: Sorted
%   - d: No equal discomforts
%   - p: p_1 > p_2 > ... > p_n
function x_inf = n_arcs_stationary_flows_kref(d,T,p,k_ref,s_max,s_bar,s_min,P_go,P_home)   
    % Parameters 
    epsl_cmp = 1e-10;
    epsl_eta = 1e-5;
    n = length(p);
    % Project p and k_ref into suitable space
    p = round(p);
    k_ref = round(k_ref);
    % Find gcd among all
    p_factor = p(1);
    for t = 1:length(p)
        p_factor = gcd(p_factor,p(t));
    end
    % Karma limits (according to Section III in [1])
    k_min = max([0 k_ref+p(end)*(T+1)]);
    k_max = k_ref+p(1)*(T+1)-p(end);
    % Karma indexing
    % k_i = k_ref + (i+i0-1)*p_factor;
    i0 = ceil((k_min-k_ref)/p_factor);
    i_max = round(((k_max-k_ref)/p_factor) + 1-i0);
    k_i = @(i) k_ref+(i+i0-1)*p_factor;
    i_p1 = find(k_i((1:i_max)') == p(1));
    % Compute gamma
    gamma = gamma_n_arcs_individual(d,T,p,k_i((1:i_max)'),k_ref,s_min,s_bar,s_max);
    % Probability of choosing j given karma k P_go*P(j|k_i)
    P_j_given_k = zeros(i_max,n);
    for j = 1:n
        if j~= 1
            P_j_given_k(:,j) = P_go*(gamma(j-1,:)-gamma(j,:))'*s_bar/s_max;
        else
            P_j_given_k(:,j) = P_go*(s_max-gamma(j,:)*s_bar)'/s_max;
        end      
    end
    % Build possibly transition matrix A
    A = P_home*eye(i_max);
    for i_1 = 1:i_max
        for i_2 = 1:i_max
            for j = 1:n
                if abs(k_i(i_1)-(k_i(i_2)-p(j))) < epsl_cmp
                    A(i_1,i_2) = A(i_1,i_2) + P_j_given_k(i_2,j);
                end
            end
        end
    end
    % Remove all unreachable components
    reachable = false(i_max,1);
    searched = false(i_max,1);
    reachable(i_p1,1) = true;
    % Search all that are reachable but not searched
    while sum(reachable.*(~searched))
        for i = 1:i_max
            if ~searched(i) && reachable(i) 
                reachable(A(:,i)~=0) = true;
                searched(i) = true;
            end
        end
    end
    A_irreducible = A(reachable,reachable);
    irreducible_idx = (1:i_max)';
    irreducible_idx = irreducible_idx(reachable);

    % Compute irreducible stationary karma distribution
    eta = zeros(i_max,1);
    eta(i_p1,1) = 1;
    eta_irreducible = eta(irreducible_idx);
    while norm((A_irreducible-eye(size(A_irreducible,1)))*eta_irreducible) > epsl_eta
        eta_irreducible = A_irreducible*eta_irreducible;
    end

    % Expand to the whole domain
    eta = zeros(i_max,1);
    eta(reachable) = eta_irreducible;
    
    % Compute stationary flows
    x_inf = (eta'*P_j_given_k)';
end

% Karma thresholds for unitary decisions, i.e., [y_bar]_j = 1 for some j
function k_ij = k_ij(i,j,k_ref,p,T)
    k_ij = k_ref+p(i)+T*p(j);
end

% Faster population generation function (CustomCreationFcn)
function population = CustomCreationFcn(GenomeLength,FitnessFcn,options,...
    x_star_pp_basis,p_norm)    
    % Init initial population matrix
    population = zeros(options.PopulationSize,GenomeLength);
    for i = 1:options.PopulationSize
        % Fill with random vectors generated by the basis perpendicular to
        % x_star
        while true         
            aux = x_star_pp_basis*(rand(GenomeLength-1,1)-0.5);
            % Make sure the prices are sorted, the first is positive, and 
            % the last negative
            if issorted(aux,'descend') && aux(1) > 0 && aux(end) < 0
                population(i,:) = round(p_norm*aux'/sqrt(aux'*aux));
                break;
            end
        end
    end
end
