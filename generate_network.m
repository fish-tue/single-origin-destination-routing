%% GENERATE_NETWORK - generate illustrative network
% Description:
%   Generates an illustrative single-destination network. Used for the
%   generation of the illustrative network in [1]. Arc ordering acording to
%   Assumption V.1 in [1].
% Outputs:
%   - network.mat (n, M, s_min, s_bar, s_max, P_home, P_go, T, alpha, beta,
%   d0, kappa, d, c0, c, x_star, d_star according to the notation in [1])
% Assumptions and limitations:
%   - Sensitivity distribution is uniform
% Other m-files required: none
% MAT-files required: none
% Toolboxes required: 
%   - CVX [2],[3]
%   - SDPT3 [4],[5]
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W. P. M. H. Heemels, and M. Salazar, “Urgency-aware 
%   optimal routing in repeated games through artificial currencies” 
%   [not published yet]
%   [2] Michael Grant and Stephen Boyd. CVX: Matlab software for 
%   disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, 
%   September 2013.
%   [3] Michael Grant and Stephen Boyd. Graph implementations for nonsmooth 
%   convex programs, Recent Advances in Learning and Control (a tribute to 
%   M. Vidyasagar), V. Blondel, S. Boyd, and H. Kimura, editors, pages 
%   95-110, Lecture Notes in Control and Information Sciences, Springer, 
%   2008.
%   [4] K.C. Toh, M.J. Todd, and R.H. Tutuncu, SDPT3 — a Matlab software 
%   package for semidefinite programming, Optimization Methods and 
%   Software, 11 (1999), pp. 545–581.
%   [5] R.H Tutuncu, K.C. Toh, and M.J. Todd, Solving semidefinite-
%   quadratic-linear programs using SDPT3, Mathematical Programming Ser. B,
%   95 (2003), pp. 189–217.

%% Initialization
clear;
rng(2); % Reproducibility

%% Parameters
% Network
n = 5; % number of links
% Population
M = 1000;
% Behaviour
% Sensitivity is an uniform distribution
s_min = 0;
s_bar = 1;
s_max = 2;
P_home = 0.05; % probability staying at home
P_go = 1 - P_home; % probability of travelling 
T = 4; % individual decision window
% Discomfort function (BPR)
alpha = 0.15;
beta = 4;
d0 = rand(n,1)*0.5+0.1;
kappa = rand(n,1)*0.5+0.1;
d = @(x) d0.*(1+alpha.*(x./kappa).^beta);
% Cost function
c0 = rand(n,1);
c = @(x) sum(c0.*d0.*(1+alpha*pow_pos(x,beta)./kappa.^beta));

%% Compute system optimum
% CVX solver
cvx_begin quiet
    % Variables
    variable x_star(n);
    % Cost funcion
    minimize(c(x_star));
    % Constraints
    ones(1,n)*x_star == P_go; % Agents that travel
    0<= x_star <= 1;
cvx_end
cost_soc_opt = cvx_optval;

%% Compute arc ordering
% d_1(x_1^star) < ... < d_n(x_n^star)
[~,srt_idx] = sort(d(x_star));
d0 = d0(srt_idx);
kappa = kappa(srt_idx);
d = @(x) d0.*(1+alpha.*(x./kappa).^beta);
c0 = c0(srt_idx);
c = @(x) sum(c0.*d0.*(1+alpha*pow_pos(x,beta)./kappa.^beta));
x_star = x_star(srt_idx);
fprintf("System opt. flows (x_star):\t");
fprintf("%g ",x_star);
fprintf("\n");
% Discomforts corresponding to socially social optimum flows
d_star = d(x_star);
fprintf("System opt. disc. (d_star):\t");
fprintf("%g ",d_star);
fprintf("\n");
% Cost corresponding to socially social optimum flows 
c_star = c(x_star);
fprintf("System opt. cost:\t\t");
fprintf("%g ",c_star);
fprintf("\n");

%% Save illustrative network
save('network.mat','n','M','s_min','s_bar','s_max','P_home','P_go','T',...
    'alpha','beta','d0','kappa','d','c0','c','x_star','d_star');
clear;