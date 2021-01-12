function complementary = complementary_new(complementary)

% Copyright (C) Harry Dankowicz, Mingwu Li
% $Id: efunc_new.m 3103 2019-06-07 17:35:49Z hdankowicz $

complementary.funcs = [];

funcarrays = { 'zero' 'embedded' 'regular' 'singular' };
for i = 1:numel(funcarrays)
  complementary.(funcarrays{i}) = [];
end

complementary.v0          = [];
complementary.tv          = [];
complementary.v_dim       =  0;
complementary.f_dim       =  0;
complementary.v_idx       = [];
complementary.f_idx       = [];
complementary.idx2par     = {};
complementary.identifyers = {};

end
