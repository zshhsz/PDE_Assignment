function [ YMat ] = Adams( func, tvec, y_init, order )
%　　Adams预测-修正算法，用于求解常微分初值问题
%   输入四个参数：函数句柄func（接收列向量、返回列向量），积分时间列向量tvec，初值行向量y_init，阶数order；
%   输出一个参数：数值解，每一行对应积分时间列向量的一行，各列为变量一个分量。
switch order
    case '4'
        row = size(tvec, 1); col = size(y_init, 2);
        YMat = zeros(row, col);
        YMat(1:4, :) = Runge_Kutta(func, tvec(1:4), y_init, '4');
        for i=4:row - 1
            stepsize = tvec(i + 1) - tvec(i);
            ydiff0 = func(tvec(i), YMat(i, :).');
            ydiff1 = func(tvec(i - 1), YMat(i - 1, :).');
            ydiff2 = func(tvec(i - 2), YMat(i - 2, :).');
            ydiff3 = func(tvec(i - 3), YMat(i - 3, :).');
            y_predict = YMat(i, :).' + (55*ydiff0 - 59*ydiff1 + 37*ydiff2 - 9*ydiff3)*stepsize/24;
            y_corrector = YMat(i, :).' + (9*func(tvec(i + 1), y_predict) + 19*ydiff0 - 5*ydiff1 + ydiff2)*stepsize/24;
            YMat(i + 1, :) = y_corrector.';
        end
end
end
