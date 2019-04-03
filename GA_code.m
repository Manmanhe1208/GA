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
for i = 1 : size(x,1)
    f(i,1) = 21.5 + x(i,1) * sin(4 * pi * x(i,1)) + x(i,2) * sin(20 * pi * x(i,2));
end
% ��100�� f ֵ�����������У�����������ֵ
[B,IX] = sort(f);
% ���100�� f ֵ�����ŵ�
paretof = B(popnum,1);
disp('�����ʼ���Ž�f��') 
disp(paretof)
% ������Ž��Ӧ��x1,x2ֵ
paretox1 = x(IX(popnum,1),1);
disp('�����ʼ����x1:')
disp(paretox1)
paretox2 = x(IX(popnum,1),2);
disp('�����ʼ����x2:')
disp(paretox2)


iter = 100;
for i = 1 : iter
    % ѡ��
    for j = 1 : popnum
        % randperm��õ��÷����Ƿ���һ����1-n�İ���n������������У�ÿ������ֻ����һ�Σ�����������������ʽ
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
    % ���100�� f ֵ�����ŵ�
    paretofp= B(popnum,1);
    if paretofp>paretof
        paretof=paretofp;
        paretox1 = parent(IX(popnum,1),1);
        paretox2 = parent(IX(popnum,1),2);
    end
     % ģ������ƽ���
     pc = 1; %�������
     hc = 20;
     offspringc = [];
     % popnum/2 ��Ϊ�˱����Ӵ���Ⱥ��ģ����100
     for j = 1 : popnum/2
	     s = randperm(popnum);
         % ����0~1 ֮��������
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
             %�������޺�����
             offsp1(1,1)=max (offsp1(1,1),xmin1);offsp1(1,1)=min (offsp1(1,1),xmax1);
             offsp1(1,2)=max (offsp1(1,2),xmin2);offsp1(1,2)=min (offsp1(1,2),xmax2);
             offsp2(1,1)=max (offsp2(1,1),xmin1);offsp2(1,1)=min (offsp2(1,1),xmax1);
             offsp2(1,2)=max (offsp2(1,2),xmin2);offsp2(1,2)=min (offsp2(1,2),xmax2);
             offspringc = [offspringc;offsp1;offsp2]; 
         else 
             %disp('������');            
             %��������ǲ������Ӵ������Ը�������
             offspringc = [offspringc;parent(s(1),:) ;parent(s(2),:)]; 
         end 
     end
     % ����
     % ����ʽ����
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
             %�������ţ�����x1,x2��ȡֵ��Χ��ͬ
             offspringm (m,1) = offspringc (m,1) + mu * (xmax1-xmin1);
              %�������޺�����
             offspringm (m,1)=max (offspringm (m,1),xmin1);
             offspringm (m,1)=min (offspringm (m,1),xmax1);
             offspringm (m,2) = offspringc (m,2) + mu * (xmax2-xmin2);
             %�������޺�����
             offspringm (m,2)=max (offspringm (m,2),xmin2);
             offspringm (m,2)=min (offspringm (m,2),xmax2);
         else
             %disp('�����б���');
             offspringm (m,1) = offspringc (m,1);
             offspringm (m,2) = offspringc (m,2);
         end
     end
     %�ϲ��Ӵ��͸�����Ⱥ
     x=[parent;offspringm];
     for j = 1 : size(x,1)
        f(j,1) = 21.5 + x(j,1) * sin(4 * pi * x(j,1)) + x(j,2) * sin(20 * pi * x(j,2));
     end   
end
disp('������Ž⣺') 
disp(paretof)
% ������Ž��Ӧ��x1,x2ֵ
paretox1 = parent(IX(popnum,1),1);
disp('�������x1:')
disp(paretox1)
paretox2 = parent(IX(popnum,1),2);
disp('�������x2:')
disp(paretox2)



