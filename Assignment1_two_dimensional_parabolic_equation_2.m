%% 用六点对称格式，ADI法，预较法和LOD法求解二维抛物方程的初边值问题
clc;
clear;
% ut = (4^-2)*(uxx + uyy)
% u(0,y,t) = u(1,y,t) = 0
% uy(x,0,t) = uy(x,1,t) = 0
% u(x,y,0) = sin(pi*x)*cos(pi*y)
%% ADI法
ax=0;ay=0;
bx=1;by=1;
N=40;
h=1/N;
x=[ax:h:bx];
y=[ay:h:by];
T=1600;
tao=1/T;
r=tao/(h^2);
a=1/16;
U=ones(N+1,N+1);
for j =1:N+1
    U(1,j)=0;
    U(N+1,j)=0;
end
for i = 2:N
    for j = 1:N+1
        U(i,j) = sin(pi*x(i))*cos(pi*y(j));
    end
end
diag_0=(1+r*a)*ones(N-1,1);
diag_1=(-r*a/2)*ones(N-2,1);
A=diag(diag_0)+diag(diag_1,1)+diag(diag_1,-1);
A2=zeros(N+1);
A2(2:N,2:N)=A;
A2(1,1)=1;
A2(N+1,N)=-1;
A2(2,1)=-r*a/2;
f=zeros(N-1,1);
f2=zeros(N+1,1);
for n=1:T
    for k=2:N
        for j=1:N-1
            f(i)=r*a/2*(U(j,k)+U(j+2,k))+(1-r*a)*U(j+1,k);
        end
        U(2:N,k)=A\f;
    end
    for j =2:N
        for k=2:N
            f2(k)=r*a/2*(U(j,k+1)+U(j,k-1))+(1-r*a)*U(j,k);
        end
        U(j,:)=(A2\f2)';
    end
end
exact=zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        exact(i,j)=sin(pi*x(i))*cos(pi*y(j))*exp(-pi^2/8);
    end
end
deta=abs(U-exact);
deta_max=max(max(deta));
fprintf('最大误差%f\n',deta_max);
figure(1);
[x_l,y_l]=meshgrid(x);
mesh(x_l,y_l,deta);
title('误差网格分布')
