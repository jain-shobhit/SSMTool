function adjoint = adjoint_new(adjoint)

adjoint.funcs  = [];
adjoint.identifyers = {};

adjoint.a_dim      = zeros(1,2);
adjoint.ax_idx     = [];
adjoint.af_idx     = [];
adjoint.l0         = [];
adjoint.tl0        = [];
adjoint.lactive    = [];
adjoint.linactive  = [];
adjoint.l_idx      = [];
adjoint.lnames     = [];
adjoint.s_idx     = [];
adjoint.snames    = [];
adjoint.sactive   = [];
adjoint.sinactive = [];
adjoint.midx       = {};

adjoint.F    = @adjoint_F;
adjoint.FDF  = @adjoint_FDF;
adjoint.DFDX = @adjoint_DFDX;

end
