clear all;close all;clc;
%% 设置参数

l = 1;
a = 1;
tmax = 0.1;
k = 0.0002;
h = 0.02;
r2 = 5/9;

%% 所求解方程
x = 0:h:l;
t = 0:k:tmax;
o = length(x);
p = length(t);
u2 = zeros(o,p);
[X.T] = meshgrid(x,t);

%u1(x,0)初始层
for i = 1:o
    if x(i) <= 1/2
        u2(i,1) = 2*x(i);
    else
        u2(i,1) = 2-2*x(i);
    end
end

%u(0,t),u(l,t)边界条件
u2(l,:) = 0;
u2(o,:) = 0;

%% 向前差分格式
for j = 1:(p-1)
    for i = 2:(o-1)
        u2(i,j+1) = r2*u2(i+1,j)+(1-2*r2)*u2(i,j)+r2*u2(i-1,j);
    end
end




