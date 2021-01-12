%% test simplified telescope formula
nmax=3;
x=sym('x',[nmax,1]);
clear ts
for n=1:nmax
    cf=(dec2bin(2^(n-1):2^n-1)-'0')*2-1;
    lfac=prod(sym(cf),2)/factorial(sym(n))/2^(sym(n)-1);
    ts(n)=simplify(lfac'*((cf*x(1:n)).^n))
end
