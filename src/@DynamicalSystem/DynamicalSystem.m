classdef DynamicalSystem < matlab.mixin.SetGetExactNames
    % DynamicalSystem construct a dynamical system object in first order or
    % second order form
    
    % M\ddot{x} + C\dot{x} + K x + fnl(x,xd) = fext(t) - Second order
    % B \dot{z} = F(z) + Fext(t)                    - First order
    % Here fnl(x) is a polynomial function of degree two or higher, which
    % is stored as a cell array such that fnl{k} corresponds to polynomials
    % of degree k+1. fnl{k} is given by a tensor of order k+2, where the
    % first mode corresponds to indices for the force vector.
    % Likewise, F(z) is a polynomial function of degree one or higher,
    % i.e., F(z) = Az + Higher order terms. F is stored as a cell array,
    % where the i-th entry gives the tensor/multiindex representation of
    % polynomials of degree i.
    
    % The second order form is converted into the first order form with z =
    % [x;\dot{x}], B = [C M;M 0], A = [-K 0;0 M], F(z)=Az+[-fnl(x,xd);0], and
    % Fext(t) = [fext(t); 0]
    
    properties
        M = []
        C = []
        K = []
        A = []
        B = []
        BinvA
        fnl = []
        F = []
        fext = []
        Fext = []
        Omega = []
        
        
        n                   % dimension for x
        N                   % dimension for z
        order = 2;          % whether second-order or first-order system
        degree              % degree of (polynomial) nonlinearity of the rhs of the dynamical system
        nKappa              % Fourier Series expansion order for Fext
        kappas =   []       % matrix with all kappas in its rows
        
        spectrum = []       % data structure constructed by linear_spectral_analysis method
        Options = DSOptions()

    end
    
    methods
        %% SET methods
        function set.A(obj,A)
            obj.A = A;
            set(obj,'order',1); % since second-order system is assumed by default
        end        

        function set.fnl(obj,fnl)
            % sets nonlinearity in second order form in multi-index format
            if iscell(fnl)
                % Input is sptensor
                % sets nonlinearity in second order form in multi-index format
                obj.fnl = set_fnl(fnl);
            else
                % Already in multi-index notation
                obj.fnl = fnl;
           end

        end
        %% GET methods
        function A = get.A(obj)
            
            if obj.order ==1
                A = obj.A;
            elseif obj.order == 2
                A = [-obj.K,         sparse(obj.n,obj.n);
                    sparse(obj.n,obj.n),   obj.M];
            end
            
        end
        
        function B = get.B(obj)
            if obj.order ==1
                
                if isempty(obj.B)
                    B = speye(obj.N,obj.N);
                else
                    B = obj.B;
                end
                
            elseif obj.order == 2
                
                B = [obj.C,    obj.M;
                    obj.M,  sparse(obj.n,obj.n)];
            end
        end
        
        function BinvA = get.BinvA(obj)
            BinvA = [sparse(obj.n,obj.n), speye(obj.n,obj.n)
                        -obj.M\obj.K,   -obj.M\obj.C];
        end
        
        function F = get.F(obj)
            
            switch obj.Options.notation
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Computation in tensor format     %%
                %   nonlinearitiy is input as tensor %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'tensor'
                    
                    if obj.order == 1
                        
                        if ~isempty(obj.F)
                            F = obj.F;
                        elseif ~isempty(obj.fnl)
                            % In this case, the DS is input to be first order, but
                            % the nonlinearity is given in second order format.
                            % See for instance example BenchmarkSSM1stOrder
                            
                            F = fnl_to_Ftens(obj);                           
                            
                        else
                            F = [];
                        end
                        
                    else
                        
                        F = fnl_to_Ftens(obj);                           

                    end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Computation in multiindex format     %%
                %   nonlinearitiy is input in multiindex %%
                %   or converted to it from tensor       %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                case 'multiindex'
                                
                    if obj.order == 1
                                                
                        if ~isempty(obj.F)
                            F = obj.F;
                        elseif ~isempty(obj.fnl)                            
                            % In this case, the DS is input to be first order, but
                            % the nonlinearity is given in second order format.
                            % See for instance example BenchmarkSSM1stOrder
                            
                            F = fnl_to_Fmulti(obj);

                        else
                            F = [];
                        end
                        
                    else  % Second order dynamical system was provided
                        
                        F = fnl_to_Fmulti(obj);
                    end
                    
                    
                otherwise
                    error('The option should be tensor or multiindex.');
            end
        end
            
        function n = get.n(obj)
            n = length(obj.M);
        end
        
        function N = get.N(obj)
            N = length(obj.A);
        end
        
        function nKappa = get.nKappa(obj)
            nKappa = numel(obj.Fext.data);
        end
        
        function kappas = get.kappas(obj)
            %kappas stored in rows
            sz_kappa = size(obj.Fext.data(1).kappa,2);
            kappas = reshape([obj.Fext.data.kappa],sz_kappa,[]).';
            
        end
        
        function Fext = get.Fext(obj)
            if obj.order ==1
                Fext = obj.Fext;
            elseif obj.order == 2
                Fext.data    = set_Fext(obj);
                Fext.epsilon = obj.fext.epsilon;
            end                           
        end
        
        
        function degree = get.degree(obj)
            degree = 0;
            if ~isempty(obj.A)
                degree = length(obj.F);
            end
        end
        
        %% other methods
        
        [V, D, W] = linear_spectral_analysis(obj)
        
        function add_forcing(obj,f,varargin)
            if ~isfield(f, 'data') % new format
                switch obj.order
                    case 1
                        nn = size(f,1);
                        Kappas = varargin{1};
                        data(1).kappa = Kappas(1);
                        data(2).kappa = Kappas(2);
                        data(1).F_n_k(1).coeffs = f(:,1);
                        data(1).F_n_k(1).ind    = zeros(1,nn);
                        data(2).F_n_k(1).coeffs = f(:,2);
                        data(2).F_n_k(1).ind    = zeros(1,nn);
                        f_ext.data = data;
                    case 2
                        nn = size(f,1);
                        Kappas = varargin{1};
                        data(1).kappa = Kappas(1);
                        data(2).kappa = Kappas(2);
                        data(1).f_n_k(1).coeffs = f(:,1);
                        data(1).f_n_k(1).ind    = zeros(1,nn);
                        data(2).f_n_k(1).coeffs = f(:,2);
                        data(2).f_n_k(1).ind    = zeros(1,nn);
                        f_ext.data = data;
                end
                
                if numel(varargin)>1
                    f_ext.epsilon = varargin{2};
                else
                    f_ext.epsilon = 1;
                end
            else
                f_ext = f;
            end
            
            switch obj.order

                case 1                    
                    obj.Fext.data = f_ext.data;                        

                    if isfield(f_ext,'epsilon')
                        obj.Fext.epsilon = f_ext.epsilon;
                        
                    elseif nargin == 3
                        obj.Fext.epsilon = varargin{1};

                    else
                        obj.Fext.epsilon = 1;

                    end

                case 2                   
                    obj.fext.data = f_ext.data;                        

                    if isfield(f_ext,'epsilon')
                        obj.fext.epsilon = f_ext.epsilon;
                        
                    elseif nargin == 3
                        obj.fext.epsilon = varargin{1};

                    else
                        obj.fext.epsilon = 1;

                    end
            end
        end

        
        fext = compute_fext(obj,t,x,xd)
        Fext = evaluate_Fext(obj,t,z)
        fnl = compute_fnl(obj,x,xd)
        dfnl = compute_dfnldx(obj,x,xd)
        dfnl = compute_dfnldxd(obj,x,xd)
        Fnl = evaluate_Fnl(obj,z)
        f = odefun(obj,t,z)
        [r, drdqdd,drdqd,drdq, c0] = residual(obj, q, qd, qdd, t)
    end
end


function [F] = fnl_to_Ftens(obj)
d = length(obj.fnl) + 1;
F = cell(1,d);
F{1} = sptensor(obj.A);

for j = 2:d
    sizej = obj.N*ones(1,j+1);
    if isempty(obj.fnl(j-1)) || isempty(obj.fnl(j-1).coeffs)
        F{j} = sptensor(sizej);
    else
        [fnl_t] = multi_index_to_tensor(obj.fnl(j-1).coeffs,obj.fnl(j-1).ind);
        subsj = fnl_t.subs;
        valsj = -fnl_t.vals;
        if obj.order==1
            valsj = -valsj;
        end
        F{j} = sptensor(subsj,valsj,sizej);
    end
    
end
end

function [F] = fnl_to_Fmulti(obj)
d = length(obj.fnl) + 1;
F = repmat(struct('coeffs',[],'ind',[]),1,d);

for j = 2:d
    if isempty(obj.fnl(j-1))
        F(j) = [];
        %% following elseif could be removed if inputs are always strictly either first
        % or second order - not like in bernoulli beam
    elseif size(obj.fnl(j-1).coeffs,1) == obj.N %fnl already 1st order form
        
        F(j).coeffs = obj.fnl(j-1).coeffs;
        F(j).ind    = obj.fnl(j-1).ind;
        
    else % conversion to 1st order form
        F(j).coeffs = [-obj.fnl(j-1).coeffs;...
            sparse(obj.n, size(obj.fnl(j-1).coeffs,2)) ];
        if obj.n == size(obj.fnl(j-1).ind,2) % No nonlinear damping
            F(j).ind = [obj.fnl(j-1).ind.';...
                sparse(obj.n, size(obj.fnl(j-1).ind,1)) ].';
        else %Nonlinear damping
            F(j).ind = obj.fnl(j-1).ind;
        end
        
    end
end
end

function [fnl_multi]  = set_fnl(fnlTensor)
%Sets second order nonlinear force in multi-index format
d   = length(fnlTensor) + 1;

fnl_multi = repmat(struct('coeffs',[],'ind',[]),1,d-1);

for j = 2:d
    if isempty(fnlTensor{j-1}) || nnz(fnlTensor{j-1}) == 0

    else
        sizej = fnlTensor{j-1}.size;
        subsj = fnlTensor{j-1}.subs;
        valsj = fnlTensor{j-1}.vals;
        tmp = tensor_to_multi_index(sptensor(subsj,valsj,sizej));
        fnl_multi(j-1).coeffs = tmp.coeffs;
        fnl_multi(j-1).ind = tmp.ind;
    end
end
end

function [data] = set_Fext(obj)
% Creates the data struct from the input second order force
% Structs for storing the coefficients
F_n_k = repmat(struct('coeffs',[],'ind',[]),numel(obj.fext.data(1).f_n_k),1);
data  = repmat(struct('kappa',[],'F_n_k',[]),numel(obj.fext.data),1);

% Fill the structs
for i = 1:numel(obj.fext.data)
    for j = 1:numel(obj.fext.data(i).f_n_k)
        F_n_k(j).coeffs = [obj.fext.data(i).f_n_k(j).coeffs;...
            sparse(obj.n, size(obj.fext.data(i).f_n_k(j).coeffs,2)) ];
        if size(obj.fext.data(i).f_n_k(j).ind,2) == obj.n
            F_n_k(j).ind = [obj.fext.data(i).f_n_k(j).ind.';...
                sparse(obj.n, size(obj.fext.data(i).f_n_k(j).ind,1)) ].';
        elseif  size(obj.fext.data(i).f_n_k(j).ind,2) == obj.N
            F_n_k(j).ind =  obj.fext.data(i).f_n_k(j).ind;
        elseif size(obj.fext.data(i).f_n_k(j).ind,2) > 0;
            error('Wrong dimensionality of external force ')
        end

    end
    data(i).F_n_k = F_n_k;
    data(i).kappa = obj.fext.data(i).kappa;
end
end