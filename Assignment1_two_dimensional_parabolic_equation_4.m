%% 用六点对称格式，ADI法，预较法和LOD法求解二维抛物方程的初边值问题
clc;
clear;
% ut = (4^-2)*(uxx + uyy)
% u(0,y,t) = u(1,y,t) = 0
% uy(x,0,t) = uy(x,1,t) = 0
% u(x,y,0) = sin(pi*x)*cos(pi*y)
%% LOD法
% main程序

function[] = LOD()
A = a^(-2);
ax = 0;bx = 1;
ay = 0;by = 1;
t0 = 1;
h = 1/40;
tao = 1/1600;
LOD_chafen(A,ax,bx,ay,by,to,h,tao)
end
%求解函数
function fT=Ture(x,y,t)
fT = sin(pi*x)*cos(pi*y)*exp(-pi^2*t/8);
end
%LOD差分函数
function[] = LOD_chafen(A,ax,bx,ay,by,to,h,tao)
tic
NX = (bx-ax)/h;
NY = (by-ay)/h;
N = NX+1;
Node N^2;
r = A*tao/(h^2);
coefM = sparse(eye(Node));
R = sparse(zeros(Node,1));
for j = 2:N-1
    for i = 2:N-1
        k = i+(j-1)*N;
        coefM(k,k-1)=-r/2;
        coefM(k.k)=1+r;
        coefM(k,k+1)=-r/2;
    end
end
Mat = sparse(zeros(Node,1));
for i = 1:N
    for j = 1:N
        Mat((i-1)*N+j)=sin(pi*(i-1)*h*cos(pi*(j-1)*h));
    end
end
for m = 1:10
    for i = 1:N
        coefM(i,i-N)=-1;
    end
    for i =2:N-1
        coefM(1+(i-1)*N,2+(i-1)*N)=0;
        coefM(i*N,i*N-1)=0;
    end
    for j = 2:N-1
        for i = 2:N-1
            R(i+(j-1)*N)=r/2*Mat(j+(i-2)*N)+r/2*Mat(j+i*N)+(1-r)*Mat(j+(i-1)*N);
        end
    end
    [Mat,z] = bicgstab(coefM,R,1e-6,100);
    
    for i=2:N-1
        for j = 1:N-1
            R(j+(i-1)*N)=r/2*Mat((j-2)*N+i)+(1-r)*Mat(i+(j-1)*N)+r/2*Mat(i+j*N);
        end
    end
    for i = 1:N
        coefM(i,i+N)=0;
    end
    for i = Node-N+1:Node
        coefM(i,i-N)=0;
    end
    for i = 2:N-1
        coefM(1+(i-1)*N,2+(i-1)*N)=-1;
        coefM(i*N,i*N-1)=-1;
    end
    
    [Mat,z]=bicgstab(coefM,R,1e-6,100);
end

lod=sparse(zeros(N,N));
true=sparse(zeros(N,N));

for i = 1:N
    for j = 1:N
        lod(i,j)=Mat(j+(i-1)*N);
        true(i,j)=True((i-1)*h,(j-1)*h,tao*10);
    end
end
for i =1:3
    for j = 1:3
        TRUE(i,j)=true(i*NX/4,j*NX/4);
        LOD(i,j)=lod(i*NX/4,j*NX/4);
    end
end
error = LOD-TRUE;
x=ax:h:bx;
y=ay:h:by;
[xx,yy]=meshgrid(x,y);
figure(1);
axis([ax bx ay by -1 1])
surf(xx,yy,(true-lod)');
title('精确解与差分解的误差');
toc
end


