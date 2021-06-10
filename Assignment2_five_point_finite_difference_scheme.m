%% 忘记要解哪一个泊松方程了，所以随便找了一个泊松方程进行五点差分格式求解
%% 先给出边界初值以及设出泊松方程所需相应变量
clc;
n=5;
h=1/n;
F=zeros((n+1)^2,1);
A=zeros(size(F));
U=ones((n+1)^2,1);  
U0=zeros((n+1)^2,1);
for i=1:n+1    
    U0(i)=0;            
    U0(n^2+n+i)=(i*h)^2;
    U0(1+(i-1)*(n+1))=0;
    U0(i*(n+1))=(i*h)^2;
end
%% 下面由五点差分格式对A进行赋值，得到五点差分法系数矩阵
for i=1:(n+1)^2
    A(i,i)=4;
end
for i=1:n-1 
    for j=2+i*(n+1):(n-1)+i*(n+1)
        A(j,j+1)=-1;
    end
end
for i=1:n-1 
    for j=3+i*(n+1):n+i*(n+1)
        A(j,j-1)=-1;
    end
end
for i=1:n-2 
    for j=2+i*(n+1):n+i*(n+1)
        A(j,j+n+1)=-1;
    end
end
for i=2:n-1 
    for j=2+i*(n+1):n+i*(n+1)
        A(j,j-n-1)=-1;
    end
end
%% 下面对F进行赋值
for i=1:n-1
    for j=2+i*(n+1):n+i*(n+1)
        F(j)=-2*((floor(j/(n+1)))^2+(mod(j,n+1)-1)^2)*h^4;
    end
end
for i=1:n-1 
    F(n+i*(n+1))=F(n+i*(n+1))+(i*h)^2;  
end
for i=1:n-1
    F((n+1)*(n-1)+i+1)=F((n+1)*(n-1)+i+1)+(i*h)^2;
end
for i=i:n+1
    F(i)=0;
end
for i=1:n-1
    F(1+i*(n+1))=0;
    F((i+1)*(n+1))=4*(i*h)^2;
end
for i=(n+1)*n+1:(n+1)^2 
    if i==(n+1)^2
        F(i)=4*(n^2)*(h^2);
    else
        F(i)=4*((mod(i,n+1)-1)^2)*(h^2);
    end
end

%% 下面对U进行赋初值
d=sum(sum(U0))/4/n;  
for i=1:(n+1)^2
    U(i)=d;
end
%% 使用高斯赛德迭代求解KU=F
eps=1e-4;  
N=500;    
u=GaussSeidel(A,F,N);  
function [U,k]=GaussSeidel(A,F,N)
U=diag(diag(A))-triu(A);
U0=zeros(length(A),1);
M=tril(A)^(-1)*U;
g=tril(A)^(-1)*F;
U=M*U0+g;
k=1;
while norm(U-U0,2)>=eps
    U0=U;
    U=M*U0+g;
    k=k+1;
    if k>=N
        break
    end
end
end

