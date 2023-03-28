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
%   - YALMIP [2]
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W. P. M. H. Heemels, and M. Salazar, “Urgency-aware 
%   optimal routing in repeated games through artificial currencies” 
%   [not published yet]
%   [2] Lofberg, Johan. "YALMIP: A toolbox for modeling and optimization in
%   MATLAB." In 2004 IEEE international conference on robotics and 
%   automation, pp. 284-289. IEEE, 2004.

%% Initialization
clear;
rng(1); % Reproducibility

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
d0 = 0.5+0.5*sort(rand(n,1));
kappa = rand(n,1);
d = @(x) d0.*(1+alpha.*(x./kappa).^beta);
% Cost function
c0 = 0.5+0.5*rand(n,1);
c = @(x) (c0.*d0.*(1+alpha.*(x./kappa).^beta))'*x;

%% Compute system optimum
% YALMIP to parse 
% Variables
x_star = sdpvar(n,1);
% Constraints
constr = [ones(1,n)*x_star == P_go;...
          0<= x_star <= 1];
% Optimize
optimize(constr,c(x_star));
% Get value
x_star = value(x_star);

%% Compute arc ordering
% d_1(x_1^star) < ... < d_n(x_n^star)
[~,srt_idx] = sort(d(x_star));
d0 = d0(srt_idx);
kappa = kappa(srt_idx);
d = @(x) d0.*(1+alpha.*(x./kappa).^beta);
c0 = c0(srt_idx);
c = @(x) (c0.*d0.*(1+alpha.*(x./kappa).^beta))'*x;
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