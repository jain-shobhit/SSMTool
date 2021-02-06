function [ r, drdqdd, drdqd, drdq] = residual_reduced_linear( q, qd, qdd, t, Assembly, Fext)
%  RESIDUAL_REDUCED_LINEAR In the following function, we construct the residual needed for time integration 
% of second-order reduced system
% 
% $\underbrace{\mathbf{V}^{\top}\mathbf{M}\mathbf{V}}_{\mathbf{M_V}}\ddot{\mathbf{q}} 
% + \underbrace{\mathbf{V}^{\top}\mathbf{C}\mathbf{V}}_{\mathbf{C_V}}\dot{\mathbf{q}} 
% + \underbrace{\mathbf{V}^{\top}\mathbf{K}\mathbf{V}}_{\mathbf{K_V}}\mathbf{q} 
% =\mathbf{V}^{\top}\mathbf{F}_{ext}(t)$,
% 
% where $\mathbf{V}\in\mathbb{R}^{n\times m}$ is the reduction basis. We use 
% the residual is defined as
% 
% $\mathbf{r}(\ddot{\mathbf{q}},\dot{\mathbf{q}},\mathbf{q}) = \mathbf{M_V}\ddot{\mathbf{q}} 
% + \mathbf{C_V}\dot{\mathbf{q}} + \mathbf{K_Vq} - \mathbf{V}^{\top}\mathbf{F}_{ext}(t)$.
% 
% A generic Residual function, whose _handle_ is passed for performing Implicit 
% Newmark and Generalized-$$\alpha$$ time integrations in this code has the following 
% syntax for linear systems:
% 
% $$\texttt{[r, drdqdd, drdqd, drdq] = Residual(q,qd,qdd,t);}$$
% 
% where
% 
% $$\texttt{r} = \mathbf{r},\\\texttt{drdqdd} = \frac{\partial\mathbf{r}}{\partial\ddot{\mathbf{q}}},\\\texttt{drdqd} 
% = \frac{\partial\mathbf{r}}{\partial\dot{\mathbf{q}}},\\\texttt{drdq}= \frac{\partial\mathbf{r}}{\partial\mathbf{q}}$$
% 
% The extra arguments:
%% 
% # $\texttt{Assembly}$, which is an instance of ReducedAssembly class
% # $\texttt{Fext}$, which is a function handle for the external forcing on 
% the full structure,
% # $\texttt{V}$ is the reduction basis $\mathbf{V}$ 
%% 
% are required for computing the residual in this case.
% 
% New residual functions that follow the above-mentioned syntax can be written 
% according to user preference. This way, the same time integration class can 
% be used to solve a variety of problems.
% 
% Please refer to the Mechanical directory in the examples folder to understand 
% applications and usage.
%% 
% In this function, it is assumed that the matrices $\mathbf{M_V,C_V,K_V}$ for 
% the finite element mesh were precomputed and stored in the $\texttt{DATA}$ property 
% of the $\texttt{Assembly}$ object to avoid unnecessary assembly during each 
% time-step. Note that this function may be modified appropriately to allow for 
% adaptive selection. i.e., when the basis changes online during time integration.
M_V = Assembly.DATA.M;
C_V = Assembly.DATA.C;
K_V = Assembly.DATA.K;
F_V = Assembly.V.' * Fext(t);
%% 
% Residual is computed according to the formula above:
r = M_V * qdd + C_V * qd + K_V * q - F_V ;
drdqdd = M_V;
drdqd = C_V;
drdq = K_V;
end