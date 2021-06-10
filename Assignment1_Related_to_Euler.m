clear all;close all; clc
f=@(x,y)1-2*x*y/(1+x^2);
h=0.5;
x=[0:h:2];
N=size(x,2)-1;
%% Exact_solver
y0=[0,zeros(0,N)];
for n=1:N
    y0(n+1)=(x(n+1)*(3+x(n+1)^2))/(3*(1+x(n+1)^2));
end
disp('精确解如下：');
disp(y0);

%% Euler method
y1=[0,zeros(1,N)];
e1=[0,zeros(1,N)];
for n=1:N
    y1(n+1)=y1(n)+h*f(x(n),y1(n));
    e1(n+1)=y1(n+1)-y0(n+1);
end
disp('欧拉方法数值解与误差如下：');
disp(y1);
disp(e1);

%% Improved Euler method
y2=[0,zeros(1,N)];
e2=[0,zeros(1,N)];
for n=1:N
    y2(n+1)=y2(n)+h*f(x(n),y2(n));
    y2(n+1)=y2(n)+h/2*(f(x(n),y2(n))+f(x(n+1),y2(n+1)));
    e2(n+1)=y2(n+1)-y0(n+1);
end
disp('改进欧拉方法数值解与误差如下：');
disp(y2);
disp(e2);

%% Classic fourth-order explicit Kutta formula 
y3=[0,zeros(1,N)];
e3=[0,zeros(1,N)];
for n=1:N
    K1=f(x(n),y3(n));
    K2=f(x(n)+1/3*h,y3(n)+1/3*h*K1);
    K3=f(x(n)+2/3*h,y3(n)-1/3*h*K1+h*K2);
    K4=f(x(n)+h,y3(n)+h*K1-h*K2+h*K3);
    y3(n+1)=y3(n)+h/8*(K1+3*K2+3*K3+K4);
    e3(n+1)=y3(n+1)-y0(n+1);
end
disp('四阶龙格-库方法数值解与误差如下：');
disp(y3);
disp(e3);

