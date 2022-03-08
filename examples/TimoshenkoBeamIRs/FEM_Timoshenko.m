function [A,K,f,x,x_dot,C,M,A_bc]=FEM_Timoshenko(n,alpha)

%This code is essentially a function which gives as output the stiffness, 
%damping and mass matrices, K, C and M, respectively, of the beam model, as
%well as the nonlinear (material response) forcing vector f, and the global
%DOF of the beam x_dot, as well as their time derivatives x_dot, and takes 
%as input the geometric and material parameters of the beam, which here are
%defined as height h, width b, length L, Young's modulus E, shear modulus 
%G, damping ratio "damping_ratio" and density rho, as well as the number of
%elements we want to divide the beam into, defined as "n". Taking out the 
%lines where these respective parameters are defined in the code below and
%deleting the "%" in the top and bottom line "transforms" this code into a
%real function. Everything below is annotated to clarify the various steps.

%Define the number of elements you want to split the beam into. 1 Element
%has 8 degrees of freedom, so entering any number "n" yields a model with
%"8*n" degrees of freedom. However, 4*(n-1) of these degrees of freedom are
%"overlapping", i.e. defined on the same nodes, so the total model will in
%the end have "8+4*(n-1)" or "4*(1+n)" degrees of freedom.


%Specify boundary conditions: 1 for "First node clamped", 2 for "Both ends
%pinned"

bc=1;

%Allocate space for stiffness, damping and mass matrices and nonlinear
%forcing vector as well as the vectors containing the global DOF's,
%x_tilde_tot and x_tilde_dot_tot

K_tot=zeros(8+4*(n-1)+n,8+4*(n-1)+n);
C_tot=zeros(8+4*(n-1)+n,8+4*(n-1)+n);
M_tot=zeros(8+4*(n-1)+n,8+4*(n-1)+n);
f_tot=zeros(8+4*(n-1)+n,1);
x_tilde_tot=f_tot;
x_tilde_dot_tot=f_tot;

%Define variable vectors "x_tilde" and its time derivative "x_tilde_dot".
%The entries of these vectors are the degrees of freedom of the beam
%element and their time derivatives.

R=[1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0 1];
J=inv(R);

%Define the beam's parameters and from those the dimensionless parameters.
%The viscoelastic constants are estimated from a damping ratio (obtained
%ideally from small amplitude vibration experiments), since this is the 
%only quantity readily available for various materials. In the current 
%example (11.4.2017), the material paramaters for steel have been used.

% SSMTOOL 1.0
% h=100;                                    %Height of beam [mm]
% b=100;                                    %Width of beam [mm]
% L=1000;                                  %Length of beam [mm]

h=40;                                    %Height of beam [mm]
b=40;                                    %Width of beam [mm]
L=1200;                                  %Length of beam [mm]
Lel=L/n;                                 %Auxiliary variable needed for integration
E=90000;                                 %Young's modulus [kPa] 
G=E/(2*1.3);                             %Shear modulus [kPa]
rho=7850*10^(-9);                        %Density [kg/mm^3]
damping_ratio=0.2;                       %Damping Ratio [-] 0.0015; 0.2 
mu=damping_ratio*2*sqrt(E*rho*h*b);      %Viscoelastic constant for axial deformation [kPa*s]
gamma=damping_ratio*2*sqrt(G*rho*h*b);   %Viscoelastic constant for shear deformation [kPa*s]

m0=b*h*rho;
m1=b*h^3/12*rho;

%Define shape functions for discretization

syms x
phi_1(x)=[1-3*x^2+2*x^3;x-2*x^2+x^3;3*x^2-2*x^3;-x^2+x^3];
phi_1=[phi_1(x/Lel)];
phi_1(x)=[1 0 0 0;0 Lel 0 0;0 0 1 0; 0 0 0 Lel]*phi_1;
phi_2(x)=[1-3*x+2*x^2;4*x-4*x^2;-x+2*x^2];
%phi_2(x)=[1-x^2;x-x^2;x^2];
phi_2(x)=[phi_2(x/Lel)];
%phi_2=[1 0 0;0 L 0;0 0 1]*phi_2;
phi_3(x)=[1-x;x];
phi_3(x)=[phi_3(x/Lel)];
phi(x)=[phi_1;phi_3;phi_2];

% In above, the local DOFs are in the form (Note that phi=[phi1;phi3;phi2]
% [u_a, u'_a, u_b, u_b', phi_a, phi_b, w_a, w_ab, w_b], where a/b/ab
% represents left/middle/right point of the element. The global DOFs is
% obtained by R*local_dofs, yielding
% [u_a, u'_a, phi_a, w_a, w_ab, u_b, u'_b, phi_b, w_b]


for i=1:n
x_tilde_full = sym('x_%d',[5*n+4,1]);
x_tilde_dot_full = sym('x_dot_',[5*n+4,1]);
x_tilde=J*x_tilde_full(5*(i-1)+1:5*i+4,:);
x_tilde_dot=J*x_tilde_dot_full(5*(i-1)+1:5*i+4,:);


%Define auxiliary vectors and matrices to enter equations more easily. A_u,
%A_w and A_phi transform the "full" DOF vectors x_tilde and x_tilde_dot
%into the smaller DOF vectors corresponding to axial, transverse and shear
%DOF.

A_w=blkdiag(zeros(4,4),zeros(2,2),eye(3,3));
A_w(1:6,:)=[];
A_phi=blkdiag(zeros(4,4),eye(2,2),zeros(3,3));
A_phi(7:9,:)=[];
A_phi(1:4,:)=[];
A_u=blkdiag(eye(4,4),zeros(5,5));
A_u(5:9,:)=[];
u_tilde=A_u*x_tilde;
u_tilde_dot=A_u*x_tilde_dot;
w_tilde=A_w*x_tilde;
w_tilde_dot=A_w*x_tilde_dot;
phi_tilde=A_phi*x_tilde;
phi_tilde_dot=A_phi*x_tilde_dot;
u=transpose(phi_1)*u_tilde;
w=transpose(phi_2)*w_tilde;
ph1=transpose(phi_3)*phi_tilde;
u_dot=transpose(phi_1)*u_tilde_dot;
w_dot=transpose(phi_2)*w_tilde_dot;
ph1_dot=transpose(phi_3)*phi_tilde_dot;

%Define entries of nonlinear material response vector
%f(x_tilde,x_tilde_dot) and stiffness matrix etc.
Mxx0lin=b*h*(E*(diff(u))+mu*diff(u_dot));
Mxx0nl=b*h*(E*(1/2*(diff(w))^2)+mu*(diff(w)*diff(w_dot)));
Mxx0=b*h*(E*(diff(u)+1/2*(diff(w))^2)+mu*(diff(u_dot)+diff(w)*diff(w_dot)));
Mxz0lin=b*h*(G*(ph1+diff(w))+gamma*(ph1_dot+diff(w_dot)));
Mxz0nl=b*h*(G*(diff(u)*ph1)+gamma*(diff(u_dot)*ph1+diff(u)*ph1_dot));
Mxz0=b*h*(G*(ph1+diff(w)+diff(u)*ph1)+gamma*(ph1_dot+diff(w_dot)+diff(u_dot)*ph1+diff(u)*ph1_dot));
Mxx1lin=b*h^3/12*(E*diff(ph1)+mu*diff(ph1_dot));
Mxx1nl=0;
Mxx1=b*h^3/12*(E*diff(ph1)+mu*diff(ph1_dot));
Mxz1lin=0;
Mxz1nl=b*h^3/12*(G*ph1*diff(ph1)+gamma*(ph1_dot*diff(ph1)+ph1*diff(ph1_dot)));
Mxz1=b*h^3/12*(G*ph1*diff(ph1)+gamma*(ph1_dot*diff(ph1)+ph1*diff(ph1_dot)));
fu=int(diff(phi_1)*(Mxx0nl+Mxz0*ph1),x,0,Lel);
fw=int(diff(phi_2)*(Mxz0nl+Mxx0*diff(w)),x,0,Lel);
fphi=int(diff(phi_3)*(Mxx1nl+Mxz1*ph1)+(phi_3)*((Mxz0nl+Mxz0*diff(u))+Mxz1*diff(ph1)),x,0,Lel);
f=[fu;fphi;fw];
M_u=m0*int(phi_1*transpose(phi_1),x,0,Lel);
M_w=m0*int(phi_2*transpose(phi_2),x,0,Lel);
M_phi=m1*int(phi_3*transpose(phi_3),x,0,Lel);
%M_3_var=m1*phi_3*transpose(phi_3);
%M_3=M_3_var(r1(1))*d1(1);
K_u=int(diff(phi_1)*(b*h*E*(diff(transpose(phi_1)))),x,0,Lel);
%K_1_var=diff(B_u*phi)*(b*h*E*(diff(transpose(B_u*phi))));
%K_1=K_1_var(r1(1))*d1(1);
C_u=mu/E*K_u;
% [r2,d2]=lgwt(2,0,Lel);
% [r1,d1]=lgwt(1,0,Lel);
%%K_2=int(diff(B_w*phi)*b*h*G*(transpose(B_phi*phi)+diff(transpose(B_w*phi))),x,0,L);
K_w=int(diff(phi_2)*b*h*G*diff(transpose(phi_2)),x,0,Lel);
K_w_phi=int(diff(phi_2)*b*h*G*(transpose(phi_3)),x,0,Lel);
%K_w_phi_var=diff(phi_2)*b*h*G*(transpose(phi_3));
%K_w_phi=K_w_phi_var(r1(1))*d1(1);

%K_2_var=diff(B_w*phi)*b*h*G*diff(transpose(B_w*phi));
%K_2_phi_var=diff(B_w*phi)*b*h*G*(transpose(B_phi*phi));
%=K_2_var(r1(1))*d1(1);
%K_2_phi=K_2_phi_var(r1(1))*d1(1);
C_w=gamma/G*K_w;
C_w_phi=gamma/G*K_w_phi;
%%K_3=int(diff(B_phi*phi)*b*h^3/12*(E*diff(transpose(B_phi*phi)))+(B_phi*phi)*b*h*G*(transpose(B_phi*phi)+diff(transpose(B_w*phi))),x,0,L);
K_phi=int(diff(phi_3)*b*h^3/12*(E*diff(transpose(phi_3)))+(phi_3)*b*h*G*(transpose(phi_3)),x,0,Lel);
K_phi_w=int((phi_3)*b*h*G*(diff(transpose(phi_2))),x,0,Lel);
%K_phi_w_var=(phi_3)*b*h*G*(diff(transpose(phi_2)));
%K_3_var=diff(phi_3)*b*h^3/12*(E*diff(transpose(phi_3)))+(phi_3)*b*h*G*(transpose(phi_3));
%K_phi=K_3_var(r1(1))*d1(1);
%K_phi_w=K_phi_w_var(r1(1))*d1(1);
%K_3_var_E=diff(B_phi*phi)*b*h^3/12*(E*diff(transpose(B_phi*phi)));
%K_3_var_G=(B_phi*phi)*b*h*G*(transpose(B_phi*phi));
%K_3_w_var=(B_phi*phi)*b*h*G*(diff(transpose(B_w*phi)));
%K_3_w=K_3_w_var(r1(1))*d1(1)+K_3_w_var(r2(2))*d2(2);
C_phi=int(diff(phi_3)*b*h^3/12*(mu*diff(transpose(phi_3)))+(phi_3)*b*h*gamma*(transpose(phi_3)),x,0,Lel);
%C_phi_w=int((phi_3)*b*h*gamma*(diff(transpose(phi_2))),x,0,L);
C_phi_w=gamma/G*K_phi_w;
%C_3_var=diff(phi_3)*b*h^3/12*(mu*diff(transpose(phi_3)))+(phi_3)*b*h*gamma*(transpose(phi_3));
%C_phi=C_3_var(r1(1))*d1(1);
%C_3_var=diff(B_phi*phi)*b*h^3/12*(mu*diff(transpose(B_phi*phi)))+(B_phi*phi)*b*h*gamma*(transpose(B_phi*phi));
%C_3_w_var=(B_phi*phi)*b*h*gamma*(diff(transpose(B_w*phi)));
%C_3=C_3_var(r1(1))*d1(1);
%C_3_w=C_3_w_var(r1(1))*d1(1);
%K_var_E(1:2,1:2)=K_1_var;
%K_var_E(3:6,3:6)=zeros(4,4);
%K_var_E(7:8,7:8)=K_3_var_E;
%K_var_G(1:2,1:2)=zeros(2,2);
%K_var_G(3:6,3:6)=K_2_var(r1(1))*d1(1);
%K_var_G(7:8,7:8)=K_3_var_G(r1(1))*d1(1);
%K_var_G(3:6,7:8)=K_2_phi_var(r1(1))*d1(1);
%K_var_G(7:8,3:6)=K_3_w_var(r1(1))*d1(1);
%K_E=K_var_E;
%K_G=K_var_G;
%K=K_E+K_G;
%C=mu/E*K_E+gamma/G*K_G;

M=blkdiag(M_u,M_phi,M_w);
K=blkdiag(K_u,K_phi,K_w);
C=blkdiag(C_u,C_phi,C_w);
K(7:9,5:6)=K_w_phi;
K(5:6,7:9)=K_phi_w;
C(7:9,5:6)=C_w_phi;
C(5:6,7:9)=C_phi_w;
%C(1:4,1:4)=C(1:4,1:4)+alpha/m0*M(1:4,1:4);
%C(7:9,7:9)=C(7:9,7:9)+alpha*M(7:9,7:9);
C=C+alpha/m0*M;


%Assemble the total model with total stiffness (K_tot), damping (C_tot) and
%mass (M_tot) matrices and global nonlinear forcing vector f_tot. The 
%remaining global degrees of freedom and their time derivatives are 
%contained in x_tilde_tot and x_tilde_dot_tot. R is a permutation matrix
%which permutes the entries of K, C, M, f, x_tilde, x_tilde_dot to 
%facilitate the assembly of the corresponding global system vectors and
%matrices. Permuted quantities are denoted by (name)_rotated. The assembly 
%occurs (as always in FEM) after the condition that the expression for the 
%virtual work of the assembled system is equal to the virtual work of all 
%the parts of the system.

x_tilde_rotated = R*x_tilde;
x_tilde_dot_rotated = R*x_tilde_dot;
f_rotated= R*f;
K_rotated= R*K*transpose(R);
C_rotated= R*C*transpose(R);
M_rotated= R*M*transpose(R);
K_tot=K_tot+blkdiag(zeros(5*(i-1),5*(i-1)),K_rotated,zeros(size(K_tot)-size(K_rotated)-size(zeros(5*(i-1),5*(i-1)))));
C_tot=C_tot+blkdiag(zeros(5*(i-1),5*(i-1)),C_rotated,zeros(size(C_tot)-size(C_rotated)-size(zeros(5*(i-1),5*(i-1)))));
M_tot=M_tot+blkdiag(zeros(5*(i-1),5*(i-1)),M_rotated,zeros(size(M_tot)-size(M_rotated)-size(zeros(5*(i-1),5*(i-1)))));
sizef=size(f_tot)-size(f_rotated);
sizef1=size(f_tot)-size(zeros(5*(i-1),1))-size(f_rotated);
sizef2=size(f_tot)-size(zeros(4+5*(i-1),1))-size(x_tilde_rotated(5:9,1));
if i==1
f_tot=f_tot+[f_rotated;zeros(sizef(1),1)];
x_tilde_tot=x_tilde_tot+[x_tilde_rotated;zeros(sizef(1),1)];
x_tilde_dot_tot=x_tilde_dot_tot+[x_tilde_dot_rotated;zeros(sizef(1),1)];
else
f_tot=f_tot+[zeros(5*(i-1),1);f_rotated;zeros(sizef1(1),1)];
x_tilde_tot=x_tilde_tot+[zeros(4+5*(i-1),1);x_tilde_rotated(5:9,1);zeros(sizef2(1),1)];
x_tilde_dot_tot=x_tilde_dot_tot+[zeros(4+5*(i-1),1);x_tilde_dot_rotated(5:9,1);zeros(sizef2(1),1)];
end
end


%Implementation of essential boundary conditions (to get rid of the
%singularities in the K and C matrices, which correspond to rigid body 
%motions): The full system Mx''+Cx'+Kx+f(x,x')=0 is transformed into a new,
%smaller reduced equivalent (up to the DOF's where boundary conditions are 
%specified) system with M_tot -> M_tot_red, C_tot -> C_tot_red, M_tot ->
%M_tot_red, f_tot -> f_tot_red, x_tilde_tot -> x_tilde_tot_red, 
%x_tilde_dot_tot -> x_tilde_dot_tot_red, where the subscript "_red" refers
%to the reduced system. A_bc is an auxiliary matrix to transform the
%system. Depending on the conditions specified, A_bc will change shape. In
%the current example (11.4.2017), the boundary conditions 

%correspond to a clamping of the end of the beam at x=0, so that all DOF's
%on the first node are set equal to zero.

if bc==1
%Clamped on first node
A_bc=[zeros(4,5*n);eye(5*n,5*n)];
A_bc=[ [0;1;0;0;zeros(5*n,1)] A_bc];
elseif bc==2
%Pinned at both ends
A_bc=zeros(4+5*n,5*n);
    if n==1
    A_bc=transpose([[0 1 0 0] zeros(1,5*n);[0 0 1 0] zeros(1,5*n);zeros(1,5*n) [0 1 0 0];zeros(1,5*n) [0 0 1 0]]);
    else 
    A_bc=A_bc+transpose([[0 1 0 0] zeros(1,5*n);[[0 0 1 0] zeros(1,5*n)];zeros(5*n-4,5*n+4);zeros(1,5*n) [0 1 0 0];zeros(1,5*n) [0 0 1 0]]);
    sizeA=size(A_bc);
    A_bc(5:sizeA(1)-4,3:sizeA(2)-2)=eye(5*n-4,5*n-4);
       % for k=1:(n-1)
       % A_bc(5+5*(k-1):4+5*k,3+5*(k-1):2+5*k)=eye(5,5);
       % end
    end
end 
  cons = [setdiff(x_tilde_tot,A_bc*transpose(A_bc)*x_tilde_tot);setdiff(x_tilde_dot_tot,A_bc*transpose(A_bc)*x_tilde_dot_tot)];
 
f_tot_red=transpose(A_bc)*subs(f_tot,cons,zeros(numel(cons),1));
x_tilde_tot_red=transpose(A_bc)*x_tilde_tot;
x_tilde_dot_tot_red=transpose(A_bc)*x_tilde_dot_tot;
M_tot_red=transpose(A_bc)*M_tot*A_bc;
K_tot_red=transpose(A_bc)*K_tot*A_bc;
C_tot_red=transpose(A_bc)*C_tot*A_bc;


%Rename everything for simplicity: K global stiffness matrix, C global
%damping matrix, M global mass matrix, f global nonlinear forcing (material
%response) vector, x global DOF, x_dot time derivatives of global DOF.
M_tot_red=vpa(M_tot_red,10);
f_tot_red=vpa(f_tot_red,10);
K_tot_red=vpa(K_tot_red,10);
C_tot_red=vpa(C_tot_red,10);

%Rename everything for simplicity: K global stiffness matrix, C global
%damping matrix, M global mass matrix, f global nonlinear forcing (material
%response) vector, x global DOF, x_dot time derivatives of global DOF.

M=double(M_tot_red);
K=K_tot_red;                                
C=C_tot_red;
f=f_tot_red;
x=x_tilde_tot_red;
x_dot=x_tilde_dot_tot_red;
A=double([zeros(5*n+1,5*n+1) eye(5*n+1,5*n+1);-M\K -M\C]);

% str='hh=matlabFunction(';
% stre=',''File'',''hh'');';
% for i=1:5*n+1
% if i==5*n+1
% str1=sprintf('f(%d)',i);
% else
% str1=sprintf('f(%d),',i);
% end
% str=strcat(str,str1);
% end
% str=strcat(str,stre);
% 
% eval(str);
end