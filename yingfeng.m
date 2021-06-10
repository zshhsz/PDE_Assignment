clc;
clear
a=1;
b=10;
L=5;
T=5;
h=0.5;
t=0.02;
r=(a*t)/h^2;
r1=(b*t)/h;
u=zeros(L/h+1,T/t+1);
time=0.02:t:T;
hang=1:L/h+1;
u(1,:)=1;
for n=1:length(time)
    for j=2:length(hang)
        if j<=10
            u(j,n+1)=r*u(j+1,n)+(1-2*r-r1)*u(j,n)+(r+r1)*u(j-1,n);
        else
            u(j,n+1)=(1-2*r-r1)*u(j,n)+(r+r1)*u(j-1,n);
        end
    end
end
[X,Y]=meshgrid(0:t:5,0:h:5);
Z=u;
mesh(X,Y,Z);
xlabel('time');
ylabel('space');
zlabel('u数值解')