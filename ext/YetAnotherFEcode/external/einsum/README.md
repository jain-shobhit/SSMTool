# Einstein summation for MATLAB
Input :
  - str : String like 'ik,kj-> ij'
  - varargin : Double called

Output :
  - Out : Double

Usage :
Matrix multiplication C = A*B
A(ik)*B(kj) = C(ij) -> C=einsum('ik,kj->ij',A,B)

Limitations :
- Works on 2016b and later
- Take only "diagonal terms" for an nd-array isn't permitted
      ex : 'iij,jk -> ik'
