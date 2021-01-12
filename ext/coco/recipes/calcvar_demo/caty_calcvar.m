function y = caty_calcvar(f, fp, p, mode)
%CATY_CALCVAR  'calcvar'-compatible encoding of derivatives of catenary Lagrangian.

switch mode
    case 'dLdf'
        y = sqrt(1+fp.^2);
    case 'dLdfp'
        y = f.*fp./sqrt(1+fp.^2);
    case 'd2Ldfdf'
        y = zeros(size(f));
    case 'd2Ldfpdf'
        y = fp./sqrt(1+fp.^2);
    case 'd2Ldfpdfp'
        y = f./(1+fp.^2).^(3/2);
end

end
