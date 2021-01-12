function [data, F] = complementary_call_F(opts, data, func, u, l, v)
% OFUNC_CALL_F       Evaluate complementary zero/monitor function
%
% See also OFUNC_CALL_DF, OFUNC_CALL_FDF

[data, F] = func.F(opts, data, u, l, v);

end
