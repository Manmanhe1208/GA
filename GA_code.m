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
for i = 1 : size(x,1)
    f(i,1) = 21.5 + x(i,1) * sin(4 * pi * x(i,1)) + x(i,2) * sin(20 * pi * x(i,2));
end
% 将100个 f 值进行升序排列，并返回索引值
[B,IX] = sort(f);
% 输出100个 f 值中最优的
paretof = B(popnum,1);
disp('输出初始最优解f：') 
disp(paretof)
% 输出最优解对应的x1,x2值
paretox1 = x(IX(popnum,1),1);
disp('输出初始最优x1:')
disp(paretox1)
paretox2 = x(IX(popnum,1),2);
disp('输出初始最优x2:')
disp(paretox2)


iter = 100;
for i = 1 : iter
    % 选择
    for j = 1 : popnum
        % randperm最常用的用法是是返回一个从1-n的包含n个数的随机排列（每个数字只出现一次）——以行向量的形式
        % s1 = randperm(100);
        s1 = randperm(size(x,1));
        se1 = s1(1);
        se2 = s1(2);
%         se1 = x(popnum-1,1)
%         se2 = x(popnum-1,2)
        if f(se1)>f(se2)
            parent(j,:) = x(se1,:);
        else
            parent(j,:) = x(se2,:);
        end 
    end
    f=[];
    for j = 1 : size(parent,1)
        f(j,1) = 21.5 + parent(j,1) * sin(4 * pi * parent(j,1)) + parent(j,2) * sin(20 * pi * parent(j,2));
    end 
    [B,IX] = sort(f);
    % 输出100个 f 值中最优的
    paretofp= B(popnum,1);
    if paretofp>paretof
        paretof=paretofp;
        paretox1 = parent(IX(popnum,1),1);
        paretox2 = parent(IX(popnum,1),2);
    end
     % 模拟二进制交叉
     pc = 1; %交叉概率
     hc = 20;
     offspringc = [];
     % popnum/2 是为了保持子代种群规模还是100
     for j = 1 : popnum/2
	     s = randperm(popnum);
         % 产生0~1 之间的随机数
         u = rand(1); 
         if u < 0.5
                 B = (2 * u)^(1/(hc + 1));
         else 
                 B = (1/(2*(1-u)))^(1/(hc + 1));
         end
         if u < pc 
             % s1 = randperm(100)
             offsp1 = 0.5 *[(1-B)* parent(s(1),:) + (1+B)* parent(s(2),:)];
             offsp2 = 0.5 *[(1-B)* parent(s(2),:) + (1+B)* parent(s(1),:)];
             %控制上限和下限
             offsp1(1,1)=max (offsp1(1,1),xmin1);offsp1(1,1)=min (offsp1(1,1),xmax1);
             offsp1(1,2)=max (offsp1(1,2),xmin2);offsp1(1,2)=min (offsp1(1,2),xmax2);
             offsp2(1,1)=max (offsp2(1,1),xmin1);offsp2(1,1)=min (offsp2(1,1),xmax1);
             offsp2(1,2)=max (offsp2(1,2),xmin2);offsp2(1,2)=min (offsp2(1,2),xmax2);
             offspringc = [offspringc;offsp1;offsp2]; 
         else 
             %disp('不交叉');            
             %不交叉就是不产生子代，所以父代保留
             offspringc = [offspringc;parent(s(1),:) ;parent(s(2),:)]; 
         end 
     end
     % 变异
     % 多项式变异
     pm = 0.5;
     hm = 20;
     offspringm = [];
     for m = 1 : popnum 
         r = rand(1);
 	     if r < 0.5
             mu = (2 * r)^(1/(hm + 1))-1;
         else
             mu = 1 - [2 * (1-r)]^(1/(hm + 1));
         end
         if r < pm
             %加索引号，由于x1,x2的取值范围不同
             offspringm (m,1) = offspringc (m,1) + mu * (xmax1-xmin1);
              %控制上限和下限
             offspringm (m,1)=max (offspringm (m,1),xmin1);
             offspringm (m,1)=min (offspringm (m,1),xmax1);
             offspringm (m,2) = offspringc (m,2) + mu * (xmax2-xmin2);
             %控制上限和下限
             offspringm (m,2)=max (offspringm (m,2),xmin2);
             offspringm (m,2)=min (offspringm (m,2),xmax2);
         else
             %disp('不进行变异');
             offspringm (m,1) = offspringc (m,1);
             offspringm (m,2) = offspringc (m,2);
         end
     end
     %合并子代和父代种群
     x=[parent;offspringm];
     for j = 1 : size(x,1)
        f(j,1) = 21.5 + x(j,1) * sin(4 * pi * x(j,1)) + x(j,2) * sin(20 * pi * x(j,2));
     end   
end
disp('输出最优解：') 
disp(paretof)
% 输出最优解对应的x1,x2值
paretox1 = parent(IX(popnum,1),1);
disp('输出最优x1:')
disp(paretox1)
paretox2 = parent(IX(popnum,1),2);
disp('输出最优x2:')
disp(paretox2)



