%========================================================================
% DESCRIPTION: 
% Finite Element model of an elastic rod with rectangular cross section
% area. Nodes are numbered from left to right, starting from one.
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
classdef FE_ElasticRod < MechanicalSystem
    properties
        L           % matrix recovering the coordinates q_full = L*q from
                    % generalized coordinates q (compatible with the
                    % constraints)
    end
    
    methods
        function obj = FE_ElasticRod(len,A,E,rho,...
                BCs,n_nodes)
            %% Constructor
            
            % Determine structural matrices
            [K,M,L] = FE_ElasticRod.structural_matrices(...
                len,A,E,rho,BCs,n_nodes);
            
            % Call parent constructor
            % NOTE THE SYSTEM HAS NO INHERENT DAMPING
            obj = obj@MechanicalSystem(M,0*M,K,[],[]);
            
            % Store matrix L
            obj.L = L;
        end
        function add_forcing(obj,inode,Fex1)
            % Determine index within generalized coordinates, associated
            % with given node
            iq = find_coordinate(obj,inode);
            
            % Add to existing force vector
            if isempty(obj.Fex1)
                obj.Fex1 = zeros(obj.n,1);
            end
            obj.Fex1(iq) = obj.Fex1(iq) + Fex1;
        end
        function add_nonlinear_attachment(obj,inode,type,varargin)
            % Determine index within generalized coordinates, associated
            % with given node and degree of freedom
            iq = find_coordinate(obj,inode);
            
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
            
            % Set force dirction
            w = zeros(obj.n,1);
            w(iq) = 1;
            obj.nonlinear_elements{i}.force_direction = w;
            
            % Make sure the description of nonlinearities is complete
            check_nonlinearities(obj);
        end
        function iq = find_coordinate(obj,inode)
            % Determine index within generalized coordinates, associated
            % with given node
            
            % Index in q_full
            iqfull = (inode-1)+1;
            
            % Determine index in q
            iq = find(obj.L(iqfull,:));
            
            % If it was not found, the DOF is not accessible (constrained)
            if isempty(iq)
                error(['Cannot apply forcing to specified node ' ...
                    numestr(inode) '.' ...
                    ' The coordinate is probably constrained.']);
            end
        end
    end
    
    methods (Static)
        function [K,M,L] = structural_matrices(len,A,...
                E,rho,BCs,n_nodes)
            %% 1D Finite Element Model of an elastic rod
            
            % Auxiliary parameters
            EA = E*A;
            rhoA = rho*A;
            
            % 1. Divide into elements and nodes
            n_elements = n_nodes-1;             % number of elements
            xn = linspace(0,len,n_nodes);       % location of nodes
            
            % 2. Specify element properties
            n_dofperelem = 2;    % number of degrees of freedom per element
            % Normalized shape functions and their derivatives
            psi = {@(xi) 1-xi, @(xi) xi};
            dpsi = {@(xi) -ones(size(xi)), @(xi) ones(size(xi))};
            % Initialize elements ('Smells like OOP spirit!')
            Elements(1:n_elements) = struct(...
                'm',zeros(n_dofperelem,n_dofperelem),...
                'k',zeros(n_dofperelem,n_dofperelem),...
                'f',zeros(n_dofperelem,1),'l',[],'le',[]);
            for e=1:n_elements
                % Association of nodes
                Elements(e).le = [e;e+1];
                
                % Element size
                Elements(e).l = xn(Elements(e).le(2)) - ...
                    xn(Elements(e).le(1));
                
                % Integration of the upper triangle of the matrices m,k
                for i=1:n_dofperelem
                    for j=1:i
                        Elements(e).m(i,j) = integral(@(xi) ...
                            rhoA*psi{i}(xi).*psi{j}(xi)*...
                            Elements(e).l,0,1);
                        Elements(e).k(i,j) = integral(@(xi) ...
                            EA*dpsi{i}(xi).*dpsi{j}(xi)/...
                            Elements(e).l,0,1);
                    end
                end
                % Supplement symmetry part
                Elements(e).m = Elements(e).m + Elements(e).m' - ...
                    diag(diag(Elements(e).m));
                Elements(e).k = Elements(e).k + Elements(e).k' - ...
                    diag(diag(Elements(e).k));
                %     For such simple elements we can also integrate analytically:
                %     Elements(e).m = rhoA*Elements(e).l/6*[2 1;1 2];
                %     Elements(e).k = EA/Elements(e).l*[1 -1;-1 1];
            end
            
            % 3. Assemble system M,K
            K = zeros(n_nodes,n_nodes); M = K;
            for e=1:n_elements
                K(Elements(e).le,Elements(e).le) = ...
                    K(Elements(e).le,Elements(e).le) + Elements(e).k;
                M(Elements(e).le,Elements(e).le) = ...
                    M(Elements(e).le,Elements(e).le) + Elements(e).m;
            end
            % Remove rows and colums in accordance with constraints B*q=0,
            % so that we take the null space, L, of B (i.e. B*L=0) as
            % reduced set of basis vectors, q = L*q_full where q_full
            % contains all DOFs of the free rod.
            L = eye(size(K,1));
            switch BCs
                case 'pinned-pinned'
                    L = L(:,2:end-1);
                case 'pinned-free'
                    L = L(:,2:end);
                case 'free-pinned'
                    L = L(:,1:end-1);
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