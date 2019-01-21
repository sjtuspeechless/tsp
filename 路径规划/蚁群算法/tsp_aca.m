
%% ����׼��
% ��ջ�������
clear all
clc

% �������м�ʱ��ʼ
t0 = clock;

%��������
points=xlsread('points_data.xlsx', 'B2:C53');

%% ������м��໥����
n = size(points,1);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(sum((points(i,:) - points(j,:)).^2));
        else
            D(i,j) = 1e-4;      %�趨�ĶԽǾ�������ֵ
        end
    end    
end

%% ��ʼ������
m = 75;                              % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                            % ����������Ҫ�̶�����
vol = 0.2;                           % ��Ϣ�ػӷ�(volatilization)����
Q = 10;                              % ��ϵ��
Heu_F = 1./D;                        % ��������(heuristic function)
Tau = ones(n,n);                     % ��Ϣ�ؾ���
Table = zeros(m,n);                  % ·����¼��
iter = 1;                            % ����������ֵ
iter_max = 100;                      % ���������� 
Route_best = zeros(iter_max,n);      % �������·��       
Length_best = zeros(iter_max,1);     % �������·���ĳ���  
Length_ave = zeros(iter_max,1);      % ����·����ƽ������  
Limit_iter = 0;                      % ��������ʱ��������

%% ����Ѱ�����·��
while iter <= iter_max
    % ��������������ϵ����λ��
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(n);
          start(i) = temp(1);
      end
      Table(:,1) = start; 
      % ������ռ�
      points_index = 1:n;
      % �������·��ѡ��
      for i = 1:m
          % �����·��ѡ��
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % �ѷ��ʵĵ㼯��(���ɱ�)
             allow_index = ~ismember(points_index,tabu);    % �μ�˵��1������ײ���
             allow = points_index(allow_index);  % �����ʵĵ㼯��
             P = allow;
             % �����֮��ת�Ƹ���
             for k = 1:length(allow)
                 P(k) = Tau(tabu(end),allow(k))^alpha * Heu_F(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % ���̶ķ�ѡ����һ�����ʵĵ�,������P���Ԫ���ܺͿ���1��ÿ����Ϣ��Ԫ�ر�ʾһ������
            Pc = cumsum(P);   
            target_index = find(Pc >= rand);    %���һ��0��1֮��ı������ҳ���Ϣ�ش�������Ԫ��
            target = allow(target_index(1));    %���̶ķ��ܱ�֤ÿ��ѭ������������и��ʽϴ��Ŀ��
            Table(i,j) = target;
         end
      end
      % ����������ϵ�·������
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);
          Limit_iter = 1; 
          
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
              Limit_iter = iter; 
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      % ������Ϣ��
      Delta_Tau = zeros(n,n);
      % ������ϼ���
      for i = 1:m
          % ��������
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-vol) * Tau + Delta_Tau;
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
end

%% �����ʾ
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
Time_Cost=etime(clock,t0);
disp(['��̾���:' num2str(Shortest_Length)]);
disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);
disp(['������������:' num2str(Limit_iter)]);
disp(['����ִ��ʱ��:' num2str(Time_Cost) '��']);

%% ��ͼ
figure(1)
plot([points(Shortest_Route,1);points(Shortest_Route(1),1)],...  %����ʡ�Է�ΪMatlab���з�
     [points(Shortest_Route,2);points(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(points,1)
    text(points(i,1),points(i,2),['   ' num2str(i)]);
end
text(points(Shortest_Route(1),1),points(Shortest_Route(1),2),'       ���');
text(points(Shortest_Route(end),1),points(Shortest_Route(end),2),'       �յ�');
xlabel('��λ�ú�����')
ylabel('��λ��������')
title(['ACA���Ż�·��(��̾���:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b')
legend('��̾���')
xlabel('��������')
ylabel('����')
title('�㷨�����켣')
