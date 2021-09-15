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
        
        spectrum = []       % data structure constructed by linear_spectral_analysis method
        Options = DSOptions()

    end
    
    methods
        %% SET methods
        function set.A(obj,A)
            obj.A = A;
            set(obj,'order',1); % since second-order system is assumed by default
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

                case 'tensor'
                    d = length(obj.fnl) + 1;
                    F = cell(1,d);
                    F{1} = sptensor(obj.A);

                    for j = 2:d
                        sizej = obj.N*ones(1,j+1);
                        if isempty(obj.fnl{j-1})
                            F{j} = sptensor(sizej);
                        else                                
                            subsj = obj.fnl{j-1}.subs;
                            valsj = -obj.fnl{j-1}.vals;
                            if obj.order==1
                                valsj = -valsj;
                            end
                            F{j} = sptensor(subsj,valsj,sizej);
                        end
                    end

                case 'multiindex'
                    d = length(obj.fnl) + 1;
                    F = cell(1,d);
                    F{1} = tensor_to_multi_index(sptensor(obj.A));

                    for j = 2:d
                        sizej = obj.N*ones(1,j+1);
                        if isempty(obj.fnl{j-1})
                            F{j} = [];
                        else                                
                            subsj = obj.fnl{j-1}.subs;
                            valsj = -obj.fnl{j-1}.vals;
                            if obj.order==1
                                valsj = -valsj;
                            end                            
                            F{j} = tensor_to_multi_index(sptensor(subsj,valsj,sizej));
                        end
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
            nKappa = size(obj.Fext.kappas,1);
        end
        
        function Fext = get.Fext(obj)
            if obj.order ==1
                Fext = obj.Fext;
            elseif obj.order == 2
                Fext.kappas = obj.fext.kappas;
                nKappas = size(Fext.kappas,1);
                Fext.coeffs = [obj.fext.coeffs;
                    sparse(obj.n, nKappas)];
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
        
        function add_forcing(obj,coeffs,kappas,varargin)
            if obj.order==1
                obj.Fext.coeffs = coeffs;
                obj.Fext.kappas = kappas;   
                if nargin == 4
                    obj.Fext.epsilon = varargin{1};
                else
                    obj.Fext.epsilon = 1;
                end                
            elseif obj.order==2
                obj.fext.coeffs = coeffs;
                obj.fext.kappas = kappas;   
                if nargin == 4
                    obj.fext.epsilon = varargin{1};
                else
                    obj.fext.epsilon = 1;
                end
            end
        end
        
        fext = compute_fext(obj,t)
        Fext = evaluate_Fext(obj,t)
        fnl = compute_fnl(obj,x,xd)
        dfnl = compute_dfnldx(obj,x,xd)
        dfnl = compute_dfnldxd(obj,x,xd)
        Fnl = evaluate_Fnl(obj,z)
        f = odefun(obj,t,z)
        [r, drdqdd,drdqd,drdq, c0] = residual(obj, q, qd, qdd, t)
    end
end


