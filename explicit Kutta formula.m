f=@(t,u)4*t*sqrt(u);
h=0.1;
y0=0;
y1=0.001;
y2=0.008;
y3=0.027;

%% Fourth-level fourth-order explicit Kutta formula
tic;
u4=[1,zeros(1,N)];
for i=1:N
    K1=f(t(i),u4(i));
    K2=f(t(i)+1/3*h,u4(i)+1/3*h*K1);
    K3=f(t(i)+2/3*h,u4(i)-1/3*h*K1+h*K2);
    K4=f(t(i)+h,u4(i)+h*K1-h*K2+h*K3);
    u4(i+1)=u4(i)+h/8*(K1+3*K2+3*K3+K4);
end
toc;
disp(['4阶4部步龙格库塔方法花费时间：',num2str(toc)]);
