%% Symbolic differentiation with coco demo bistable - Generation of right-hand side
% This file generates the right-hand side and (by default) its first two
% derivatives for the periodically forced Duffing oscillator with hardening
% nonlinearity (see coco demo bistable).  See also <demo.html> for
% follow-up demo. Consider the Duffing oscillator
%
% $$
%   \frac{\mathrm{d}}{\mathrm{d} t}\vec{x}=\vec{f}(t,\vec{x},\vec{p})
% $$ 
%
% where 
%
% $$
%  \vec{x}=
%   \left[\begin{array}{c}
%     x\\ v
%   \end{array}\right],\quad 
%   \vec{p}=\left[
%   \begin{array}{c}
%     T\\ a\\ \gamma
%   \end{array}\right],\quad
%   \vec{f}=\left[
%   \begin{array}{l}
%     v\\ -\gamma v-x-x^3+a \cos(2 \pi t/T)
%   \end{array}\right]
% $$
%% Load path and package (if octave is used)
% The generation of right-hand sides is octave compatible. If octave is
% used, one may need to load the package |symbolic|
clear
addpath([pwd(),'/../../toolbox']); % path of symcoco routines
if sco_isoctave()
    pkg load symbolic   % if octave is used load package symbolic
end
%% Create symbols for time, state, parameters
% Below are the standard way of declaring symbols and defining a
% symbolic expression |f| using Matlab's symbolic toolbox.
syms t x v gam a T 
f=[v; -gam*v-x-x^3+a*cos(2*pi*t/T)]; 
%% Generate code: output is side effect, written  to file
% The call to function |sco_sym2funcs| below takes the symbolic expression
% |f| of size |2x1| as the argument that will be converted into a function.
% It computes (by default) its first two derivatives and saves it into a
% Matlab function that can be called later by |coco|. The second and third
% input of |sco_sym2funcs| define the generated function as a function of 3
% arguments. The second input defines the dimension and symbols of these
% inputs, the third assigns names (in this case |'t'|, |'x'| and |'p'| such
% that this call generates a function |f(t,x,p)|, where |x| is the |2x1|
% vector |[x;v]| from the symbolic expression and |p| is the |3x1| vector
% |[T;a;gam]| (such that |p(1)=T, p(2)=a, p(3)=gam|). All other inputs of
% |sco_sym2funcs| are optional name-value pairs. Optional input
% |'filename'| determines which file the result is written into.
sco_sym2funcs(f,...             % symbolic expression for f
    {t,[x;v],[T;a;gam]},...     % which symbols are in which inputs of f
    {'t','x','p'},...           % names for inputs of f
    'vector',[0,1,1],...        % are inputs scalar or vectors (default: vectors)
    'filename','sym_bistable'); % filename for result
rmpath([pwd(),'/../../toolbox']); % path of symcoco routines
%% Demo for usage of generated functions
% Go to <demo.html> for follow-up demo.