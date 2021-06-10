%%Adams内插法
function [X,Y]=Adams4nei(fun,x0,b,y0,h)
x=x0;
y=y0;
p=128;
n=fix((b-x0)/h);
if n<5
return
end
X=zeros(p,1);
Y=zeros(p,length(y));
f=zeros(p,1);
k=1;
X(k)=x;
Y(k,:)=y';
%RK4求初值
for k=2:4
b2=1/2;b3=1/2;b4=1;
x1=x+1/2*h;x2=x+1/2*h;x3=x+h;k1=feval(fun,x,y);
y1=y+b2*h*k1;k2=feval(fun,x1,y1);
y2=y+b3*h*k2;k3=feval(fun,x2,y2);
y3=y+b4*h*k3;k4=feval(fun,x3,y3);
y=y+1/6*h*(k1+2*k2+2*k3+k4);
x=x+h;
X(k)=x;
Y(k,:)=y;
k=k+1;
end
X;Y;f(1:4)=feval(fun,X(1:4),Y(1:4));
%内插公式
for k=4:n
f(k+1)=feval(fun,X(k),Y(k));
X(k+1)=X(1)+h*k;
Y(k+1)=Y(k)+(h/24)*((f(k-2:k+1))'*[1 -5 19 9]');
f(k+1)=feval(fun,X(k+1),Y(k+1));
f(k)=f(k+1);
k=k+1;
end
X=X(1:n+1);Y=Y(1:n+1);n=1:n+1;
%%Adams外插法
function [X,Y]=Adams4wai(fun,x0,b,y0,h)
x=x0;
y=y0;
p=128;
n=fix((b-x0)/h);
if n<5
return
end
X=zeros(p,1);
Y=zeros(p,length(y));
f=zeros(p,1);
k=1;
X(k)=x;
Y(k,:)=y';
%RK4求初值
for k=2:4
b2=1/2;b3=1/2;b4=1;
x1=x+1/2*h;x2=x+1/2*h;x3=x+h;k1=feval(fun,x,y);
y1=y+b2*h*k1;k2=feval(fun,x1,y1);
y2=y+b3*h*k2;k3=feval(fun,x2,y2);
y3=y+b4*h*k3;k4=feval(fun,x3,y3);
y=y+1/6*h*(k1+2*k2+2*k3+k4);
x=x+h;
X(k)=x;
Y(k,:)=y;
k=k+1;
end
X;Y;f(1:4)=feval(fun,X(1:4),Y(1:4));
%外插公式
for k=4:n
f(k)=feval(fun,X(k),Y(k));
X(k+1)=X(1)+h*k;
Y(k+1)=Y(k)+(h/24)*((f(k-3:k))'*[-9 37 -59 55]');
f(k+1)=feval(fun,X(k+1),Y(k+1));
f(k)=f(k+1);
k=k+1;
end
X=X(1:n+1);Y=Y(1:n+1);n=1:n+1;
%%Euler法
function [X,Y]=Euler2(fun,x0,b,y0,h)
x=x0;
n=fix((b-x)/h);
X=zeros(n+1,1);
y=y0;
Y=zeros(n+1,1);
k=1;
X(k)=x;
Y(k)=y';
for k=2:n+1
X(k)=x+(k-1)*h;
fxy=feval(fun,x,y);
Y(k)=y+h*fxy;
y=Y(k);
k=k+1;
end
%%dy/dx=f(x,y)型微分方程
function dy=fun(x,y)
dy=-5*y;
%%四种方法结算结果与解析解比较，fun为f(x,y),x0，y0满足初值条件y(x0)=y0，b是x的上界，h为步长
function MyFun(fun,x0,b,y0,h)
[x1,y1]=Adams4wai(fun,x0,b,y0,h);
[x2,y2]=Adams4nei(fun,x0,b,y0,h);
[x3,y3]=Euler2(fun,x0,b,y0,h);
[x4,y4]=Euler2gaijin(fun,x0,b,y0,h);
y=dsolve('Dy=-5*y','y(0)=1','t');%解析解
plot(x1,y1,'r*','markersize',10)
hold on
plot(x2,y2,'r.','markersize',10)
hold on
plot(x3,y3,'o','markersize',10)
hold on
plot(x4,y4,'h','markersize',10)
hold on
ezplot(y,[0 1])
hold on
title('三阶Adams内插法、外插法、Euler法及改进Euler法解初值问题结果比较')
legend('Adams外插法','Adams内插法','Euler法','改进Euler法','精确解')
