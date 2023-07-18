%========================================================================
% DESCRIPTION: 
% Finite Element model of an Euler-Bernoulli beam.
% Nodes are numbered from left to right, starting from one.
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.1 Copyright (C) 2019  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
classdef FE_EulerBernoulliBeam < MechanicalSystem
    % Elastic beam oriented horizontally
    properties
        L           % matrix recovering the coordinates q_full = L*q from
                    % generalized coordinates q (compatible with the
                    % constraints)
    end
    
    methods
        function obj = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
                BCs,n_nodes)
            %% Constructor
            
            % Determine structural matrices
            [K,M,L] = FE_EulerBernoulliBeam.structural_matrices(...
                len,height,thickness,E,rho,BCs,n_nodes);
            
            % Call parent constructor
            % NOTE THE SYSTEM HAS NO INHERENT DAMPING
            obj = obj@MechanicalSystem(M,0*M,K,[],[]);
            
            % Store matrix L
            obj.L = L;
        end
        function add_forcing(obj,inode,dof,Fex1)
            % Determine index within generalized coordinates, associated
            % with given node and degree of freedom
            iq = find_coordinate(obj,inode,dof);
            
            % Add to existing force vector
            if isempty(obj.Fex1)
                obj.Fex1 = zeros(obj.n,1);
            end
            obj.Fex1(iq) = obj.Fex1(iq) + Fex1;
        end
        function add_nonlinear_attachment(obj,inode,dof,type,varargin)
            % Determine index within generalized coordinates, associated
            % with given node and degree of freedom
            [iq,dof] = find_coordinate(obj,inode,dof);
            
            % Define nonlinearity
            if isempty(obj.nonlinear_elements)
                i = 1;
            else
                i = length(obj.nonlinear_elements)+1;
            end
            obj.nonlinear_elements{i} = struct;
            for k=1:length(varargin)/2
                obj.nonlinear_elements{i}.(varargin{2*k-1}) = ...
                    varargin{2*k};
            end
            obj.nonlinear_elements{i}.type = type;
            if strcmp(dof,'rot')
                error(['No nonlinearities defined for rotational' ...
                    'degrees of freedom yet.']);
            else
                % Set force dirction
                w = zeros(obj.n,1);
                w(iq) = 1;
                obj.nonlinear_elements{i}.force_direction = w;
            end
            
            % Make sure the description of nonlinearities is complete
            check_nonlinearities(obj);
        end
        function [iq,dof] = find_coordinate(obj,inode,dof)
            % Determine index within generalized coordinates, associated
            % with given node and degree of freedom
            switch lower(dof)
                case {'translation','translational','trans','displacement'}
                    idof = 1;
                    dof = 'trans';
                case {'rotation','rotational','rot','slope'}
                    idof = 2;
                    dof = 'rot';
                otherwise
                    error(['Unexpected degree of freedom specifier ' ...
                        dof '.']);
            end
            iqfull = 2*(inode-1)+idof;  % index in q_full
            
            % Determine index in q
            iq = find(obj.L(iqfull,:));
            
            % If it was not found, the DOF is not accessible (constrained)
            if isempty(iq)
                error(['Cannot apply forcing to specified node ' ...
                    numestr(inode) ' to DOF ' dof '.' ...
                    ' The coordinate is probably constrained.']);
            end
        end
    end
    
    methods (Static)
        function [K,M,L] = structural_matrices(length,height,thickness,...
                E,rho,BCs,n_nodes)
            %% 1D Finite Element Model of an Euler-Bernoulli beam
            
            % Auxiliary parameters
            EI = E*thickness*height^3/12;
            rhoA = rho*thickness*height;
            
            % 1. Divide into elements and nodes
            Ne = n_nodes-1;                         % number of elements
            xn = linspace(0,length,n_nodes);        % location of nodes
            
            % 2. Specify element properties
            Nde = 4;
            % Normalized shape functions and their derivatives
            psi = {@(xi,l) 1-3*xi.^2+2*xi.^3, ...
                @(xi,l) l*xi.*(1-xi).^2, ....
                @(xi,l) xi.^2.*(3-2*xi), ...
                @(xi,l) l*xi.^2.*(xi-1)};
%             dpsi = {@(xi,l) -6*xi+6*xi.^2, ...
%                 @(xi,l) l*(1-xi).^2 - 2*l*xi.(1-xi), ....
%                 @(xi,l) 2*xi.*(3-2*xi) - 2*xi.^2, ...
%                 @(xi,l) l*2*xi.*(xi-1) + l*xi.^2};
            ddpsi = {@(xi,l) -6*ones(size(xi))+12*xi, ...
                @(xi,l) l*(-4*ones(size(xi))+6*xi), ....
                @(xi,l) 6*ones(size(xi))-12*xi, ...
                @(xi,l) l*(-2*ones(size(xi))+6*xi)};
            % Initialize elements ('Smells like OOP spirit!')
            Elements(1:Ne) = struct('m',zeros(Nde,Nde),'k',zeros(Nde,Nde),...
                'f',zeros(Nde,1),'l',[],'le',[]);
            for e=1:Ne
                % Association of nodes
                Elements(e).le = [e;e+1];
                
                % Element size
                len_e = xn(Elements(e).le(2)) - xn(Elements(e).le(1));
                Elements(e).l = len_e;
                
                % Integration of the upper triangle of the matrices m,k
                for i=1:Nde
                    for j=1:i
                        Elements(e).m(i,j) = integral(@(xi) ...
                            rhoA*psi{i}(xi,len_e).*psi{j}(xi,len_e)*...
                            len_e,0,1);
                        Elements(e).k(i,j) = integral(@(xi) ...
                            EI*ddpsi{i}(xi,len_e).*ddpsi{j}(xi,len_e)/...
                            len_e^3,0,1);
                    end
                    Elements(e).f(i) = integral(@(xi) psi{j}(xi,len_e)*...
                        len_e,0,1);
                end
                % Supplement symmetric part
                Elements(e).m = Elements(e).m + Elements(e).m' - ...
                    diag(diag(Elements(e).m));
                Elements(e).k = Elements(e).k + Elements(e).k' - ...
                    diag(diag(Elements(e).k));
                % For such simple elements we can also integrate analytically:
                %     Elements(e).m = rhoA*len_e/420*[...
                %         156 22*len_e 54 -13*len_e;...
                %         22*len_e 4*len_e^2 13*len_e -3*len_e^2;...
                %         54 13*len_e 156 -22*len_e;...
                %         -13*len_e -3*len_e^2 -22*len_e 4*len_e^2];
                %     Elements(e).k = EI/len_e^3*[...
                %         12 6*len_e -12 6*len_e;...
                %         6*len_e 4*len_e^2 -6*len_e 2*len_e^2;...
                %         -12 -6*len_e 12 -6*len_e;...
                %         6*len_e 2*len_e^2 -6*len_e 4*len_e^2];
            end
            
            % 3. Assemble system M,K
            K = zeros(Nde/2*n_nodes,Nde/2*n_nodes); M = K;
            for e=1:Ne
                ILE = Nde/2*(kron(Elements(e).le(:),ones(Nde/2,1))-1) + ...
                    repmat((1:Nde/2)',2,1);
                K(ILE,ILE) = K(ILE,ILE) + Elements(e).k;
                M(ILE,ILE) = M(ILE,ILE) + Elements(e).m;
            end
            
            % Remove rows and colums in accordance with constraints B*q=0,
            % so that we take the null space, L, of B (i.e. B*L=0) as
            % reduced set of basis vectors, q = L*q_full where q_full
            % contains all DOFs of the free beam.
            L = eye(size(K,1));
            switch BCs
                case 'clamped-clamped'
                    L = L(:,3:end-2);
                case 'clamped-pinned'
                    L = L(:,[3:end-2 end]);
                case 'clamped-free'
                    L = L(:,3:end);
                case 'pinned-clamped'
                    L = L(:,2:end-2);
                case 'pinned-pinned'
                    L = L(:,[2:end-2 end]);
                case 'pinned-free'
                    L = L(:,2:end);
                case 'free-clamped'
                    L = L(:,1:end-2);
                case 'free-pinned'
                    L = L(:,[1:end-2 end]);
                case 'free-free'
                    % Nothing to constrain
                otherwise
                    error('Not implemented yet.');
            end
            % Apply constraints
            K = L'*K*L;
            M = L'*M*L;
        end
    end
end