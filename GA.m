% 随机初始化

% 生成两列的上下限相等
% popnum = 100;
% xmin = -3;
% xmax = 12.1;
% for i = 1 : popnum      % 控制种群规模
%        for j = 1 : 2    % 控制有几个x
%            x(i,j) = xmin + rand * (xmax-xmin)
%        end
% end
clear
clc


popnum = 100;
xmin1 = -3;
xmax1 = 12.1;
xmin2 = 4.1;
xmax2 = 5.8;
for i = 1 : popnum
        x(i,1) = xmin1 + rand * (xmax1 -xmin1);
        x(i,2) = xmin2 + rand * (xmax2-xmin2);
end

% 输出当前种群中的最优x1,x2,以及f
for i = 1 : popnum 
    f(i,1) = 21.5 + x(i,1) * sin(4 * pi * x(i,1)) + x(i,2) * sin(20 * pi * x(i,2));
end
% 将100个 f 值进行升序排列，并返回索引值
[B,IX] = sort(f);
% 输出100个 f 值中最优的
paretof = B(popnum,1);
disp('输出当前最优解：') 
disp(paretof)
% 输出最优解对应的x1,x2值
paretox1 = x(IX(popnum,1),1);
disp('输出当前最优x1:')
disp(paretox1)
paretox2 = x(IX(popnum,1),2);
disp('输出当前最优x2:')
disp(paretox2)


iter = 50;
for i = 1 : iter
    % 二元锦标赛选择
    for j = 1 : popnum
        % randperm最常用的用法是是返回一个从1-n的包含n个数的随机排列（每个数字只出现一次）——以行向量的形式
        % s1 = randperm(100);
        s1 = randperm(size(x,1));
        se1 = s1(1);
        se2 = s1(2);
        %         se1 = x(popnum-1,1)
        %         se2 = x(popnum-1,2)
        if se1 >se2
            parent(j,:) = x(se1,:);
        else
            parent(j,:) = x(se2,:);
        end
    end
    f = [];
    for j = 1:popnum
        f(j,1) = 21.5 + parent(j,1)*sin(4*pi*parent(j,1)) + parent(j,2)*sin(20*pi*parent(j,2));
    end
    [B,IX] = sort(f);
    % 输出100个 f 值中最优的
    paretofp = B(popnum,1);
    if paretofp > paretof
        paretof = paretofp;
        paretox1 = parent(IX(popnum,1),1);
        paretox2 = parent(IX(popnum,1),2);
    end
    % 模拟二进制交叉
    pc = 0.8; %交叉概率
    hc = 20;
    offspringc = [];
    % popnum/2 是为了保持子代种群规模还是100
    for j = 1 : popnum/2
        %这里要重新产生1-n的包含n个数的随机排列，这里的n是100，上面的s1从第二次迭代开始是1-200，所以要区分，不能直接用前面的
        s = randperm(popnum);
        se1 = s(1);
        se2 = s(2);
        % 产生0~1 之间的随机数
        u = rand(1);
        if u < 0
            B = (2 * u)^(1/(hc + 1));
        else
            B = (1/(2*(1-u)))^(1/(hc + 1));
        end
        if u < pc
            % s1 = randperm(100)
            % 第一个子代交叉                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        s1 = randperm(size(x,1)
            offsp1 = 0.5 *[(1-B)* parent(se1,:) + (1+B)* parent(se2,:)];
            % 控制第一个个体上下限
            % 比如一个个体的范围是[4.1,5.8]，交叉或变异之后个体是2.1，小于下限，所以把这个个体直接替换为下限4.1，
            % 还有如果交叉或变异之后个体是8，大于上限，所以把这个个体直接替换为上限5.8。
            % offsp1(1,1)范围是[xmin1,xmax1],offsp1(1,2)范围是[xmin2,xmax2]
            offsp1(1,1) = max(offsp1(1,1),xmin1);
            offsp1(1,1) = min(offsp1(1,1),xmax1);
            offsp1(1,2) = max(offsp1(1,2),xmin2);
            offsp1(1,2) = min(offsp1(1,2),xmax2);
            % 第二个子代交叉
            offsp2 = 0.5 *[(1-B)* parent(se2,:) + (1+B)* parent(se1,:)];
            % 控制第二个个体的上下限
            % offsp2(1,1)范围是[xmin1,xmax1],offsp2(1,2)范围是[xmin2,xmax2]
            offsp2(1,1) = max(offsp2(1,1),xmin1);
            offsp2(1,1) = min(offsp2(1,1),xmax1);
            offsp2(1,2) = max(offsp2(1,2),xmin2);
            offsp2(1,2) = min(offsp2(1,2),xmax2);
            offspringc = [offspringc;offsp1;offsp2];
        else
            % 不交叉，即父代
            offspringc = [offspringc;parent(se1,:);parent(se2,:)];
        end
    end
    % 变异
    % 多项式变异
    pm = 0.5;
    hm = 20;
    offspringm = [];
    for j = 1 : popnum
        % rand 生成[0,1]区间的均匀分布随机数或随机矩阵           
        r = rand;
        if r < 0.5
            mu = (2 * r)^(1/(hm + 1))-1;
        else
            mu = 1 - [2 * (1-r)]^(1/(hm + 1));
        end
        if r < pm
            %加索引号，由于x1,x2的取值范围不同
            offspringm1= offspringc(j,1) + mu * (xmax1-xmin1);
            % 标注上下限
            offspringm1 = max(offspringc(j,1),xmin1);
            offspringm1 = min(offspringc(j,1),xmax1);
            % 标注上下限
            offspringm2=  offspringc(j,2) + mu * (xmax2-xmin2);
            offspringm2 = max(offspringc(j,2),xmin2);
            offspringm2 = min(offspringc(j,2),xmax2);
            %上面的offspringm1和offspringm2都对应一个实数，如果进行下面的操作相当于把100*2的矩阵转变为了200*1的矩阵，所以不能直接这样操作
            offspringm = [offspringm;[offspringm1,offspringm2]];
        else
            % 不进行变异
            %上面的offspringm1和offspringm2都对应一个实数，如果进行下面的操作相当于把100*2的矩阵转变为了200*1的矩阵，所以不能直接这样操作
            offspringm = [offspringm;[offspringc(j,1),offspringc(j,2)]];
        end
    end
    % 合并父代和子代种群
    x = [parent;offspringm];
    for j = 1 : size(x,1)
        f(j,1) = 21.5 + x(1,1)*sin(4*pi*x(1,1)) + x(1,2)*sin(20*pi*x(1,2));
    end
end
%输出最优解是在整个迭代循环之后的，因为我们在循环中记录最优解了，最后循环结束输出最优解就可以了
disp('输出当前最优f值:')
disp(paretof)
% 输出当前最优x1,x2
paretox1 = parent(IX(popnum,1),1);
disp('输出当前最优x1:')
disp(paretox1)
paretox2 = parent(IX(popnum,1),2);
disp('输出当前最优x2:')
disp(paretox2)



