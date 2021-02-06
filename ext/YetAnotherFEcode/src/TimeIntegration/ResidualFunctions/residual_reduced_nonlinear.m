function [ r, drdqdd,drdqd,drdq, c0] = residual_reduced_nonlinear( q, qd, qdd, t, Assembly, Fext)
%  RESIDUAL_REDUCED_NONLINEAR In the following function, we construct the residual needed for time integration 
% of
% 
% second-order system
% 
% $\mathbf{M_V}\ddot{\mathbf{q}} + \mathbf{C_V}\dot{\mathbf{q}} + \mathbf{F_V}(\mathbf{q}) 
% =\mathbf{V}^{\top}\mathbf{F}_{ext}(t)$,
% 
% where $\mathbf{M_V}:=\mathbf{V}^{\top}\mathbf{M}\mathbf{V} ,~\mathbf{C_V}:=\mathbf{V}^{\top}\mathbf{C}\mathbf{V} 
% ,~\mathbf{F_V}(\mathbf{q}):=\mathbf{V}^{\top}\mathbf{F}(\mathbf{Vq}) $ are reduced 
% operators. We use the reduced residual is defined as
% 
% $\mathbf{r}(\ddot{\mathbf{q}},\dot{\mathbf{q}},\mathbf{q}) = \mathbf{M_V}\ddot{\mathbf{q}} 
% + \mathbf{C_V}\dot{\mathbf{q}} + \mathbf{F_V}(\mathbf{q}) - \mathbf{V}^{\top}\mathbf{F}_{ext}(t)$.
% 
% A generic Residual function, whose _handle_ is passed for performing Implicit 
% Newmark and Generalized- $$\alpha$$  nonlinear time integration schemes in this 
% code has the following syntax:
% 
% $$\texttt{[r, drdqdd, drdqd, drdq, c0] = Residual(q,qd,qdd,t);}$$
% 
% where
% 
% $$\texttt{r} = \mathbf{r},\\\texttt{drdqdd} = \frac{\partial\mathbf{r}}{\partial\ddot{\mathbf{q}}},\\\texttt{drdqd} 
% = \frac{\partial\mathbf{r}}{\partial\dot{\mathbf{q}}},\\\texttt{drdq}= \frac{\partial\mathbf{r}}{\partial\mathbf{q}},\\\texttt{c0} 
% = \textrm{a scalar measure for comparing the residual norm while checking for 
% convergence}$$
% 
% The extra arguments:
%% 
% # $\texttt{Assembly}$, which is an instance of ReducedAssembly class
% # $\texttt{Fext}$, which is a function handle for the external forcing,
% # $\texttt{V}$ is the basis used for Galerkin projection $\mathbf{V}$ and 
% is a property of the ReducedAssembly Class
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
% In this function, it is assumed that the matrices $\mathbf{M_V,C_V}$ for the 
% finite element mesh were precomputed and stored in the $\texttt{DATA}$ property 
% of the $\texttt{Assembly}$ object to avoid unnecessary assembly during each 
% time-step.
V = Assembly.V;
M_V = Assembly.DATA.M;
C_V = Assembly.DATA.C;
%% 
% The reduced tangent stiffness $\mathbf{K_V} = \frac{\partial\mathbf{F_V}}{\partial\mathbf{q}}$ 
% and the the reduced internal force $\mathbf{F_V}$, however, need to be assmbled 
% at each iteration depending on the current state. This assembly is best perfomed 
% at an element level in its reduced form as
% 
% $\mathbf{K_V}(\mathbf{q}) = \sum_{e=1}^{n_e} \mathbf{V}_e^{\top}\mathbf{K}_e(\mathbf{V}_e\mathbf{q})\mathbf{V}_e,\\\mathbf{F_V} 
% (\mathbf{q}) = \sum_{e=1}^{n_e} \mathbf{V}_e^{\top}\mathbf{F}_e(\mathbf{V}_e\mathbf{q})$. 
% 
% We first obtain the full degrees of freedom vector $\mathbf{u}$ over the mesh 
% from the reduced variables $\mathbf{q}$ as $\mathbf{u} = \mathbf{Vq}$, and then 
% directly assemble the reduced nonlinear operators  with the Assembly class method 
% $\texttt{tangent\_stiffness\_and\_force\_modal(u,V)}$. 
u = V*q;
[K_V, F_V] = Assembly.tangent_stiffness_and_force(u);
%% 
% Residual is computed according to the formula above:
F_inertial = M_V * qdd;
F_damping = C_V * qd;
F_ext_V =  V.'*Fext(t);
r = F_inertial + F_damping + F_V - F_ext_V ;
drdqdd = M_V;
drdqd = C_V;
drdq = K_V;
%% 
% We use the following measure to comapre the norm of the residual $\mathbf{r}$
% 
% $$\texttt{c0} = \|\mathbf{M_V}\ddot{\mathbf{q}}\| + \|\mathbf{C_V}\dot{\mathbf{q}}\| 
% + \|\mathbf{F_V}(\mathbf{q})\| + \|\mathbf{V}^{\top}\mathbf{F}_{ext}(t)\|$$
c0 = norm(F_inertial) + norm(F_damping) + norm(F_V) + norm(F_ext_V);
end