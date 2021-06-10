%% 用六点对称格式，ADI法，预较法和LOD法求解二维抛物方程的初边值问题
clc;
clear;
format long;
% ut = (4^-2)*(uxx + uyy)
% u(0,y,t) = u(1,y,t) = 0
% uy(x,0,t) = uy(x,1,t) = 0
% u(x,y,0) = sin(pi*x)*cos(pi*y)
%% 预较法
J=40;
N=1600;
h=1/J;
t=1/N;
r=1;
a=1/16;
[U]=zeros(J+1,J+1,N+1);
[U1]=zeros(J+1,J+1,N+1);
for n=1:N+1
    for i = 1:J+1
        for j= 1:J+1
            U1(i,j,n)=sin(pi*(i-1)*h)*cos(pi*(j-1)*h)*exp(-pi^2*(n-1)*t/8);
        end
    end
end

for j=1:J+1
    for k=1:J+1
        U(j,k,1)=sin(pi*((j-1)*h))*cos(pi*((k-1)*h));
    end
end
I=ones(1,J+1);
I=I*(-a*r/2);
v=I;
u=ones(1,J+1);
for i = 2:J
    u(1,i)=1+a*r;
end
b=zeros(1,J+1);
b1=zeros(1,J+1);
y=zeros(1,J+1);
x=zeros(1,J+1);
y1=zeros(1,J+1);
x1=zeros(1,J+1);
u(1,1)=u(1,1);
for i = 2:J+1
    I(1,i)=I(1,i)/u(1,i-1);
    u(1,i)=u(1,i)-I(1,i)*v(1,i-1);
end

for n =2:N+1
    u1=zeros((J+1)*(J+1),1);
    for k=1:J+1
        b(1,i)=U(i,k,n-1);
    end
    y(1,1)=b(1,1);
    for i=2:J+1
        y(1,i)=b(1,i)-I(1,i)*y(1,i-1);
    end
    x(1,J+1)=y(1,J+1)/u(1,J+1);
    u1((J+1)*k,1)=x(1,J+1);
    for i=J:-1:1
        x(1,i)=(y(1,i)-v(1,i)*x(1,i+1))/u(1,i);
        u1((J+1)*k-(J+1-i),1)=x(1,i);
    end
end

u2=zeros((J+1)*(J+1),1);
for k=1:J+1
    for i =1:J+1
        b1(1,i)=u1(k+(J+1)*(i-1),1);
    end
    y1(1,1)=b1(1,1);
    for i=2:J+1
        y1(1,i)=b1(1,i)-I(1,i)*y1(1,i-1);
    end
    x1(1,J+1)=y1(1,J+1)/u(1,J+1);
    u2((J+1)*k,1)=x1(1,J+1);
    for i =J:-1:1
        x1(1,i)=y1(1,J+1)/u(1,J+1);
        u2((J+1)*k,1)=x1(1,J+1);
        for i=J:-1:1
            x1(1,i)=(y1(1,i)-v(1,i)*x1(1,i+1))/u(1,i);
            u2((J+1)*k-(J+1-i),1)=x1(1,i);
        end
    end
    
    for i=1:J+1
        U(1,i,n)=0;
        U(J+1,i,n)=0;
    end
    for j=3:J
        for k=3:J
            U(j,k,n)=U(j,k,n-1)+r*a(u2(k+(J+1)*j,1)+u2(k+(J+1)*(j-2),1)+u2(k+1+(J+1)*(j-1),1)+u2(k-1+(J+1)*(j-1),1)-4*u2(k+(J+1)*(j-1),1));
        end
        U(j,1,n)=U(j,2,n);
        U(j,J+1,n)=U(j,J,n);
    end
end
for i=1:3
    UTRUE(i,j)=Ut(i*10+1,j*10+1);
    PrU(i,j)=UU(i*10+1,j*10+1);
end
Errors=PrU-UTRUE;
format long;
UTRUE';
PrU';
format short;
Errors