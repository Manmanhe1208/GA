% �����ʼ��

% �������е����������
% popnum = 100;
% xmin = -3;
% xmax = 12.1;
% for i = 1 : popnum      % ������Ⱥ��ģ
%        for j = 1 : 2    % �����м���x
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

% �����ǰ��Ⱥ�е�����x1,x2,�Լ�f
for i = 1 : popnum 
    f(i,1) = 21.5 + x(i,1) * sin(4 * pi * x(i,1)) + x(i,2) * sin(20 * pi * x(i,2));
end
% ��100�� f ֵ�����������У�����������ֵ
[B,IX] = sort(f);
% ���100�� f ֵ�����ŵ�
paretof = B(popnum,1);
disp('�����ǰ���Ž⣺') 
disp(paretof)
% ������Ž��Ӧ��x1,x2ֵ
paretox1 = x(IX(popnum,1),1);
disp('�����ǰ����x1:')
disp(paretox1)
paretox2 = x(IX(popnum,1),2);
disp('�����ǰ����x2:')
disp(paretox2)


iter = 50;
for i = 1 : iter
    % ��Ԫ������ѡ��
    for j = 1 : popnum
        % randperm��õ��÷����Ƿ���һ����1-n�İ���n������������У�ÿ������ֻ����һ�Σ�����������������ʽ
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
    % ���100�� f ֵ�����ŵ�
    paretofp = B(popnum,1);
    if paretofp > paretof
        paretof = paretofp;
        paretox1 = parent(IX(popnum,1),1);
        paretox2 = parent(IX(popnum,1),2);
    end
    % ģ������ƽ���
    pc = 0.8; %�������
    hc = 20;
    offspringc = [];
    % popnum/2 ��Ϊ�˱����Ӵ���Ⱥ��ģ����100
    for j = 1 : popnum/2
        %����Ҫ���²���1-n�İ���n������������У������n��100�������s1�ӵڶ��ε�����ʼ��1-200������Ҫ���֣�����ֱ����ǰ���
        s = randperm(popnum);
        se1 = s(1);
        se2 = s(2);
        % ����0~1 ֮��������
        u = rand(1);
        if u < 0
            B = (2 * u)^(1/(hc + 1));
        else
            B = (1/(2*(1-u)))^(1/(hc + 1));
        end
        if u < pc
            % s1 = randperm(100)
            % ��һ���Ӵ�����                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        s1 = randperm(size(x,1)
            offsp1 = 0.5 *[(1-B)* parent(se1,:) + (1+B)* parent(se2,:)];
            % ���Ƶ�һ������������
            % ����һ������ķ�Χ��[4.1,5.8]����������֮�������2.1��С�����ޣ����԰��������ֱ���滻Ϊ����4.1��
            % ���������������֮�������8���������ޣ����԰��������ֱ���滻Ϊ����5.8��
            % offsp1(1,1)��Χ��[xmin1,xmax1],offsp1(1,2)��Χ��[xmin2,xmax2]
            offsp1(1,1) = max(offsp1(1,1),xmin1);
            offsp1(1,1) = min(offsp1(1,1),xmax1);
            offsp1(1,2) = max(offsp1(1,2),xmin2);
            offsp1(1,2) = min(offsp1(1,2),xmax2);
            % �ڶ����Ӵ�����
            offsp2 = 0.5 *[(1-B)* parent(se2,:) + (1+B)* parent(se1,:)];
            % ���Ƶڶ��������������
            % offsp2(1,1)��Χ��[xmin1,xmax1],offsp2(1,2)��Χ��[xmin2,xmax2]
            offsp2(1,1) = max(offsp2(1,1),xmin1);
            offsp2(1,1) = min(offsp2(1,1),xmax1);
            offsp2(1,2) = max(offsp2(1,2),xmin2);
            offsp2(1,2) = min(offsp2(1,2),xmax2);
            offspringc = [offspringc;offsp1;offsp2];
        else
            % �����棬������
            offspringc = [offspringc;parent(se1,:);parent(se2,:)];
        end
    end
    % ����
    % ����ʽ����
    pm = 0.5;
    hm = 20;
    offspringm = [];
    for j = 1 : popnum
        % rand ����[0,1]����ľ��ȷֲ���������������           
        r = rand;
        if r < 0.5
            mu = (2 * r)^(1/(hm + 1))-1;
        else
            mu = 1 - [2 * (1-r)]^(1/(hm + 1));
        end
        if r < pm
            %�������ţ�����x1,x2��ȡֵ��Χ��ͬ
            offspringm1= offspringc(j,1) + mu * (xmax1-xmin1);
            % ��ע������
            offspringm1 = max(offspringc(j,1),xmin1);
            offspringm1 = min(offspringc(j,1),xmax1);
            % ��ע������
            offspringm2=  offspringc(j,2) + mu * (xmax2-xmin2);
            offspringm2 = max(offspringc(j,2),xmin2);
            offspringm2 = min(offspringc(j,2),xmax2);
            %�����offspringm1��offspringm2����Ӧһ��ʵ���������������Ĳ����൱�ڰ�100*2�ľ���ת��Ϊ��200*1�ľ������Բ���ֱ����������
            offspringm = [offspringm;[offspringm1,offspringm2]];
        else
            % �����б���
            %�����offspringm1��offspringm2����Ӧһ��ʵ���������������Ĳ����൱�ڰ�100*2�ľ���ת��Ϊ��200*1�ľ������Բ���ֱ����������
            offspringm = [offspringm;[offspringc(j,1),offspringc(j,2)]];
        end
    end
    % �ϲ��������Ӵ���Ⱥ
    x = [parent;offspringm];
    for j = 1 : size(x,1)
        f(j,1) = 21.5 + x(1,1)*sin(4*pi*x(1,1)) + x(1,2)*sin(20*pi*x(1,2));
    end
end
%������Ž�������������ѭ��֮��ģ���Ϊ������ѭ���м�¼���Ž��ˣ����ѭ������������Ž�Ϳ�����
disp('�����ǰ����fֵ:')
disp(paretof)
% �����ǰ����x1,x2
paretox1 = parent(IX(popnum,1),1);
disp('�����ǰ����x1:')
disp(paretox1)
paretox2 = parent(IX(popnum,1),2);
disp('�����ǰ����x2:')
disp(paretox2)



