%% 用六点对称格式，ADI法，预较法和LOD法求解二维抛物方程的初边值问题
clc;
clear;
% ut = (4^-2)*(uxx + uyy)
% u(0,y,t) = u(1,y,t) = 0
% uy(x,0,t) = uy(x,1,t) = 0
% u(x,y,0) = sin(pi*x)*cos(pi*y)
%% 六点对称格式
% main程序
a1 = 0;
b1 = 1;
m = 40;
n = 40;
th = 1600;
tao = 1/th;
a = 4;
n1 = input('输入要计算的时间层n1=');
[u] = six(a,a1,b1,m,n,th,n1);
% 计算精确解
h = (b1-a1)/m;
for j = 1:m+1
    for k=1:n+1
        uu(j,k) = uexact((j-1)*h,(k-1)*h,n1*tao);
    end
end
%在节点（xj,yk）= （j/4,k/4),j,k=1,2,3的计算结果
for j = 1:3
    for k = 1:3
        Exact(j,k) = uu(j*m/4,k*n/4);%精确解
        Difference_solution(j,k) = u(j*m/4,k*n/4);%差分解
    end
end

%作图
x = a1:h:b1;
y = a1:h:b1;
[xx,yy]=meshgrid(x,y);
%画出精确解图像
figure(1)
surf(xx,yy,uu')
title('精确解')
%画出差分图像
figure(2)
surf(xx,yy,u')
title('差分解')
%画出精确解与差分解之间的误差图
figure(3)
surf(xx,yy,(uu-u)')
title('精确解与差分解的误差')
function [u]=six(a,a1,b1,m,n,th,n1)
m = 40;
n = 40;
tao = 1/th;
knots = (m+1)*(n+1);
h = 1/m;
r = tao/(a*a*h*h);
u0=[];
for j = 1:m+1
    for k = 1:n+1
        u0(j,k) = uexact((j-1)*h,(k-1)*h,0);
    end
end
%形成系数矩阵和右端项
XS = sparse(eye(knots));
Rhs = sparse(zeros(knots,1));
solution = sparse(zeros(knots,1));
%边界
for j = 1:n+1
    l = j;
    XS(l,l) = 1;
    Rhs(l) = 0;
end
for j = 1:n+1
    l=m*(n+1)+j;
    XS(l,l) = 1;
    Rhs(l) = 0;
end
for j = 1:m+1
    l = (j-1)*(n+1)+1;
    XS(l,l) = 1;
    XS(l,l+1) = -1;
end
for j = 1:m+1
    l = j*(n+1);
    XS(l,l) = 1;
    XS(l,l-1) = -1;
end
%内点编号
for j = 2:m
    for k = 2:n
        l = (j-1)*(n+1)+k;
        XS(l,l) = 1+2*r;
        XS(l,l-1) = -r/2;
        XS(l,l+1) = -r/2;
        XS(l,l-(n+1)) = -r/2;
        XS(l,l+(n+1)) = -r/2;
    end
end
%开始迭代计算
for t = 1:n1
    for j = 2:m
        for k = 2:n
            l = (j-1)*(n+1)+k;
            Rhs(l) = r/2*(u0(j+1,k)+u0(j,k+1)+u0(j-1,k)+u0(j,k-1)+(1-2*r)*u0(j,k));
        end
    end
    solution = bicgstab(XS,Rhs,1.0e-6,400);
    for j = 1:m+1
        for k = 1:n+1
            l = (j-1)*(n+1)+k;
            u(j,k) = solution(l);
        end
    end
    u0 = u;
end
end
%精确解
function [f]=uexact(x,y,t)
f = sin(pi*x)*cos(pi*y)*exp(-pi*pi*t/8);
end
