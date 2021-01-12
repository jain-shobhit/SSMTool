%% Test symbolic differentiation for coco, demo sphere_optim
% See |coco_folder/core/examples/sphere_optim| and CORE-tutorial for original |coco|
% demo, and <demo.html> for outputs of demo produced with the derivatives
% generated below.
%
% This demo shows in <demo.html> how one can use symbolic derivatives of up
% to second order for coco computations if the expected function format is
% 
%   function [data,y_or_J]=func(prob,data,u)
%
%%
clear
addpath([pwd(),'/../../toolbox']);
if sco_isoctave()
    pkg load symbolic
end
%% Introduce symbolic variables for constraint and objectuve functional
p=sym('p',[3,1]); % p for constraint
syms x            % x for constraint
%% Define symbolic expressions for constraint (sphere) and objective
constraint = x^2+p(1)^2+p(2)^2+p(3)^2-1;
obj=sum([x;p]);
%% generate code for constraint
sco_sym2funcs(...
    constraint,...
    {x,p},...
    {'x','p'},...
    'filename','sym_sphere_constraint');
%% generate code for objective
sco_sym2funcs(...
    obj,...
    {[x;p]},...
    {'u'},...
    'filename','sym_sphere_obj');
%% Remove path to |symcoco|
rmpath([pwd(),'/../../toolbox']);
