function [a,b,delta,delta1,a_ijkl,b_ijkl,c_ijkl,g_ijkl,intphi] = cal_parameters(N,l,BC)

syms x
% N = 2; % number of modes
% l = 1;

% construct modal functions
phis   = [];
dphis  = [];
ddphis = [];
dddphis  = [];
ddddphis = [];
lamda = zeros(N,1);
if strcmp(BC,'simply-simply') 
    lamda(1:N)=(1:N)*pi;%simply-supported 这个公�?得到的特�?值是精准解
elseif strcmp(BC,'clamped-clamped') 
     lamda(1)=4.730040744862704;lamda(2)=7.853204624095837;
     lamda(3)=10.995607838001671;lamda(4)=14.137165491257464;
     lamda(5)=17.278759657399480;lamda(6)=20.420352245626059;
     lamda(7)=23.561944902040455;lamda(8)=26.703537555508188;
     lamda(9)=29.845130209103253;lamda(10)=32.986722862692822;
     lamda(11)=36.128315516282626;lamda(12)=39.269908169872416;
%  lamda(2:N)=(2:N)*pi+0.5*pi;%clamped-clamped
%  这个公�?得到的特�?值是近似解，精度�?高，直接求特�?方程的解比较好
elseif strcmp(BC,'clamped-free') 
    lamda(1)=1.8751040687119611664453082410782141625701117335311;
    lamda(2)=4.6940911329741745764363917780198120493898967375458;
    lamda(3)=7.8547574382376125648610085827645704578485419292300;
    lamda(4)=10.995540734875466990667349107854702939612972774652;
    lamda(5)=14.137168391046470580917046812551772068603076792975;
    lamda(6)=17.278759532088236333543928414375822085934519635550;
    lamda(7)=20.420352251041250994415811947947837046137288894544;
    lamda(8)=23.561944901806443501520253240198075517031265990051;
    lamda(9)=26.703537555518298805445478738088174009227496691444;
    lamda(10)=29.845130209102817263788730907430405063369916050237;
    lamda(11)=32.986722862692838446168314183853683000367036041646;
    lamda(12)=36.128315516282621834281162204648446323089539352502;
% lamda(3:N)=(3:N)*pi-0.5*pi;% 这个公�?得到的特�?值是近似解，精度�?高，直接求特�?方程的解比较好
end
intphi = zeros(N,1);
for n=1:N
   if strcmp(BC,'simply-simply') 
   phin = sqrt(2)*sin(lamda(n)*x); %simply-supported
   elseif strcmp(BC,'clamped-clamped') 
   phin = cos(lamda(n)*x)-cosh(lamda(n)*x)-(cos(lamda(n))-cosh(lamda(n)))/(sin(lamda(n))-sinh(lamda(n)))*...
         (sin(lamda(n)*x)-sinh(lamda(n)*x));%clamped-clamped
   elseif strcmp(BC,'clamped-free') 
   phin = cos(lamda(n)*x)-cosh(lamda(n)*x)-(cos(lamda(n))+cosh(lamda(n)))/(sin(lamda(n))+sinh(lamda(n)))*...
        (sin(lamda(n)*x)-sinh(lamda(n)*x));%clamped-free
   end
   dphin  = diff(phin,x);
   ddphin = diff(dphin,x);
   dddphin = diff(ddphin,x);
   ddddphin = diff(dddphin,x);
   phis   = [phis; phin];
   dphis  = [dphis; dphin];
   ddphis = [ddphis; ddphin];
   dddphis = [dddphis; dddphin];
   ddddphis = [ddddphis; ddddphin];
   normphi = int(phin^2, x, 0, l);
   intphi(n) = double(int(phin, x, 0, l));
   fprintf('L2 norm of phi%d is %d\n', n, normphi);

end

% compute linear coeffficients
delta=zeros(N);delta1=zeros(N);a=zeros(N);b=zeros(N);
% for i=1:N
%     for j=1:N
%         delta(i,j) = int(phis(i)*phis(j),x,0,l);
%         delta1(i,j) = int(phis(i)*ddddphis(j),x,0,l);
%         a(i,j) = int(phis(i)*dphis(j),x,0,l);
%         b(i,j) = int(phis(i)*ddphis(j),x,0,l);
%     end
% end

if strcmp(BC,'simply-simply') 
% compute nonlinear coeffficients of two ends supported pipe
% alpha1 = zeros(N);
% alpha2 = zeros(N);
% alpha  = zeros(N,N,N,N);
% 
% for i=1:N
%     for j=1:N
%         alpha1(i,j) = int(phis(i)*ddphis(j),x,0,l);
%         alpha2(i,j) = int(dphis(i)*dphis(j),x,0,l);
%     end
% end
%         
% for n=1:N
%     for m=1:N
%         for p=1:N
%             for q=1:N
%                 alpha(n,m,p,q) = alpha1(n,q)*alpha2(m,p);
%             end
%         end
%     end
% end

filename=strcat('pinned pinned coeff_Matrix\NN=',num2str(N));
n1=load(strcat(filename,'\n1.txt'));
alpha  = zeros(N,N,N,N);
k=0;
for n=1:N
    for m=1:N
        k=k+1;
        alpha(:,:,n,m) = n1(:,1+(k-1)*N:k*N);
    end
end

% compute nonlinear coeffficients of cantilevered pipe
% 悬臂�?的特�?函数比较�?�?�，在求int（f，x,0,x）的时候�?么计算很慢，�?么无法算出。
%因此需�?改�?��?路，�?用matlab进行积分�?作，而是在maple中计算出�?�，然�?�在matlab中读�?�

elseif strcmp(BC,'clamped-free') 
    
filename=strcat('coeff_Matrix\NN=',num2str(N));
n1=load(strcat(filename,'\n1.txt'));n2=load(strcat(filename,'\n2.txt'));
n3=load(strcat(filename,'\n3.txt'));n4=load(strcat(filename,'\n4.txt'));
n5=load(strcat(filename,'\n5.txt'));n6=load(strcat(filename,'\n6.txt'));  
n7=load(strcat(filename,'\n7.txt'));n8=load(strcat(filename,'\n8.txt'));
n9=load(strcat(filename,'\n9.txt'));n10=load(strcat(filename,'\n10.txt'));
n11=load(strcat(filename,'\n11.txt'));n12=load(strcat(filename,'\n12.txt'));
n13=load(strcat(filename,'\n13.txt'));n14=load(strcat(filename,'\n14.txt'));
n15=load(strcat(filename,'\n15.txt'));n16=load(strcat(filename,'\n16.txt'));  
n17=load(strcat(filename,'\n17.txt'));n18=load(strcat(filename,'\n18.txt'));
n19=load(strcat(filename,'\n19.txt'));n20=load(strcat(filename,'\n20.txt'));
n21=load(strcat(filename,'\n21.txt'));n22=load(strcat(filename,'\n22.txt'));
n23=load(strcat(filename,'\n23.txt'));n24=load(strcat(filename,'\n24.txt'));
n25=load(strcat(filename,'\n25.txt'));
delta=load(strcat(filename,'\delta.txt'));delta1=load(strcat(filename,'\delta1.txt'));
a=load(strcat(filename,'\a1.txt'));b=load(strcat(filename,'\b1.txt'));
N1 = zeros(N,N,N,N);N2 = zeros(N,N,N,N);N3 = zeros(N,N,N,N);N4 = zeros(N,N,N,N);N5 = zeros(N,N,N,N);
N6 = zeros(N,N,N,N);N7 = zeros(N,N,N,N);N8 = zeros(N,N,N,N);N9 = zeros(N,N,N,N);N10 = zeros(N,N,N,N);
N11 = zeros(N,N,N,N);N12 = zeros(N,N,N,N);N13 = zeros(N,N,N,N);N14 = zeros(N,N,N,N);N15 = zeros(N,N,N,N);
N16 = zeros(N,N,N,N);N17 = zeros(N,N,N,N);N18 = zeros(N,N,N,N);N19 = zeros(N,N,N,N);N20 = zeros(N,N,N,N);
N21 = zeros(N,N,N,N);N22 = zeros(N,N,N,N);N23 = zeros(N,N,N,N);N24 = zeros(N,N,N,N);N25 = zeros(N,N,N,N);
% a_ijkl= zeros(N,N,N,N);c_ijkl = zeros(N,N,N,N);b_ijkl = zeros(N,N,N,N);
for n=1:N
    for m=1:N
        s=(n-1)*N+m;
        dn=N^2;
        e=s+(N-1)*dn;
        N1(:,:,n,m)=n1(:,s:dn:e);N2(:,:,n,m)=n2(:,s:dn:e);N3(:,:,n,m)=n3(:,s:dn:e);
        N4(:,:,n,m)=n4(:,s:dn:e);N5(:,:,n,m)=n5(:,s:dn:e);N6(:,:,n,m)=n6(:,s:dn:e);
        N7(:,:,n,m)=n7(:,s:dn:e);N8(:,:,n,m)=n8(:,s:dn:e);N9(:,:,n,m)=n9(:,s:dn:e);
        N10(:,:,n,m)=n10(:,s:dn:e);N11(:,:,n,m)=n11(:,s:dn:e);N12(:,:,n,m)=n12(:,s:dn:e);
        N13(:,:,n,m)=n13(:,s:dn:e);N14(:,:,n,m)=n14(:,s:dn:e);N15(:,:,n,m)=n15(:,s:dn:e);
        N16(:,:,n,m)=n16(:,s:dn:e);N17(:,:,n,m)=n17(:,s:dn:e);N18(:,:,n,m)=n18(:,s:dn:e);
        N19(:,:,n,m)=n19(:,s:dn:e);N20(:,:,n,m)=n20(:,s:dn:e);N21(:,:,n,m)=n21(:,s:dn:e);
        N22(:,:,n,m)=n22(:,s:dn:e);N23(:,:,n,m)=n23(:,s:dn:e);N24(:,:,n,m)=n24(:,s:dn:e);
        N25(:,:,n,m)=n25(:,s:dn:e);
    end
end
a_ijkl=3*N4+N5+N13-N17-N25;
c_ijkl=N2-N9-N15+N22;
b_ijkl=N1-N7-N14+N19;
g_ijkl=N6-N18;
end
end
    
                
                