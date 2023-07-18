% nlvib_decode
%
% Syntax:
% s = nlvib_decode(X, Solinfo, Sol, analysis, method, n, H)
%
% Description: this function reorganizes the raw results coming from
% solve_and_continue.m (from NLvib) in a more user-friendly format.
% INPUTS
%   X: solution from "solve_and_continue"
%   Solinfo: second output from "solve_and_continue". If not needed, just
%   put Solinfo=[].
%   Sol: third output from "solve_and_continue", which is the output of the
%   optional postprocessing function. If you're doing the shooting method
%   and you are using "nlvib_ShootingPostprocess" as postprocessing
%   function with solve_and_continue, then the fields of Sol are merged
%   (otherwise you would have Sol(i).fields for every frequency...). If not
%   needed, just put Sol=[].
%   analysis: 'FRF' or 'NMA'
%   method: 'HB' or 'SH'
%   n: number of DOFs
%   H (optional): number of harmonics. Not required for method='SH'
% OUTPUTS
%   s: a struct variable with results stored in fields with (hopefully)
%   intuitive names. The fields change depending on the method/analysis.
%
% Author: Jacopo Marconi
% Created: 26 April 2021
% Last modified: 27 April 2021

function s = nlvib_decode(X, Solinfo, Sol, analysis, method, n, H)

if nargin < 7
    H = 0;
    if strcmpi(method,'HB')
        error('Select H>0.')
    end
end

CASE = upper([analysis '-' method]);
s = struct();

switch method
    case 'HB'
        switch analysis
            case 'NMA'
                % Interpret solver output
                Psi = X(1:end-3,:);	% Normalized solution. 
                                    % Each column contains: 2H+1 
                                    % coefficients for the ndofs, i.e.:
                                    % Psi(    1 :   n, j) --> A0
                                    % Psi(  n+1 : 2*n, j) --> A1 (Re)
                                    % Psi(2*n+1 : 3*n, j) --> B1 (Im)
                                    % ...
                                    % Psi((2*H-1)*n+1 : 2*H*n, j) --> AH
                                    % Psi(2*H*n+1 : (2*H+1)*n, j) --> BH
                                    % for the j-th frequency
                
                om  = X(end-2,:);  	% circular frequency
                del = X(end-1,:);  	% modal damping ratio
                log10a = X(end,:); 	% log10 of the reference a(imod)
                a = 10.^log10a;
                
                % de-normalize solution
                Q = Psi .*repmat(a, size(Psi, 1), 1);
            case 'FRF'
                om = X(end,:);
                Q = X(1:end-1,:);
        end
        
        % reorganize Q, separating Re/Im by harmonic number
        nw = length(om);
        Q0 = Q(1:n, :);
        Qre = zeros(n, nw, H);
        Qim = zeros(n, nw, H);
        for h = 1 : H
            Qre(:, :, h) = Q((2*h-1)*n+1 : 2*h*n, :);
            Qim(:, :, h) = Q((2*h)*n+1 : (2*h+1)*n, :);
        end
        
        % store results in a struct
        s.omega = om;
        if strcmpi(analysis,'NMA')
            s.damp = del;
            s.a = a;
        end
        s.Q0 = Q0;
        s.Qre = Qre;
        s.Qim = Qim;
        s.case = CASE;
        s.legend = 'Qre/Qim: [nDOFs * nOmega * H]';
        s.Sol = Sol;
        
    case 'SH'
        switch analysis
            case 'NMA'
        
                om  = X(end-2,:);                   % rad/s
                del = X(end-1,:);                   % modal damping ratio
                log10qsinorm_shoot = X(end,:);
                a = 10.^log10qsinorm_shoot;

                s.omega = om;
                s.delta = del;
                s.a = a;
                s.X = X;        
                s.Xinfo   = ['X = [q0-/a, u0-/(omega*a), omega, delta, log10(a)]. '...
                    'q0- (u0-) is the INITIAL CONDITION for displacement (velocity) for all BUT imod-th coordinate (-)'];
                s.case = CASE;
        
        case 'FRF'
        
            om  = X(end,:);                     % rad/s

            s.omega = om;
            s.X = X;        
            s.Xinfo   = ['X = [q0, u0/omega, omega]. '...
                'q0 (u0) is the INITIAL CONDITION for displacement (velocity)'];
            s.case = CASE;
        end
        
        if ~isempty(Sol)
            try
                % reorganize Sol content for ease of use
                % WORKS ONLY if Sol has been obtained from the 
                % "nlvib_ShootingPostprocess" function!
                for ii = 1 : length(Sol)
                    so = Sol(ii);
                    s.harmonics.Qabs.displ(:,ii,:) = so.harmonics.Qabs.displ;
                    s.harmonics.Qabs.vel(:,ii,:)   = so.harmonics.Qabs.vel;
                    s.harmonics.Qphase.displ(:,ii,:) = so.harmonics.Qphase.displ;
                    s.harmonics.Qphase.vel(:,ii,:)   = so.harmonics.Qphase.vel;
                    s.harmonics.rms.displ(:, ii) = so.harmonics.rms.displ;
                    s.harmonics.rms.vel(:, ii)   = so.harmonics.rms.vel;
                    s.harmonics.max.displ(:, ii) = so.harmonics.rms.displ;
                    s.harmonics.max.vel(:, ii)   = so.harmonics.max.vel;
                    s.time.displ(:, ii, :) = so.time.displ;
                    s.time.vel(:, ii, :)   = so.time.vel;
                    s.floquet.mucrit(ii) = so.floquet.mucrit;
                    s.floquet.stable(ii) = so.floquet.stable;
                    s.floquet.ctime(ii) = so.floquet.ctime;
                end
                s.harmonics.info = so.harmonics.info;
                s.time.Nsamples = so.time.Nsamples;
                s.time.info = so.time.info;
            catch
                % ... otherwise it simply returns Sol
                s.Sol = Sol;
            end
        end
        
    otherwise
        error('Wrong analysis/method selected.')
end

s.Solinfo = Solinfo;

