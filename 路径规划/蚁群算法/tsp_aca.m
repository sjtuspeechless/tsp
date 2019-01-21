
%% 数据准备
% 清空环境变量
clear all
clc

% 程序运行计时开始
t0 = clock;

%导入数据
points=xlsread('points_data.xlsx', 'B2:C53');

%% 计算城市间相互距离
n = size(points,1);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(sum((points(i,:) - points(j,:)).^2));
        else
            D(i,j) = 1e-4;      %设定的对角矩阵修正值
        end
    end    
end

%% 初始化参数
m = 75;                              % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
vol = 0.2;                           % 信息素挥发(volatilization)因子
Q = 10;                              % 常系数
Heu_F = 1./D;                        % 启发函数(heuristic function)
Tau = ones(n,n);                     % 信息素矩阵
Table = zeros(m,n);                  % 路径记录表
iter = 1;                            % 迭代次数初值
iter_max = 100;                      % 最大迭代次数 
Route_best = zeros(iter_max,n);      % 各代最佳路径       
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度  
Limit_iter = 0;                      % 程序收敛时迭代次数

%% 迭代寻找最佳路径
while iter <= iter_max
    % 随机产生各个蚂蚁的起点位置
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(n);
          start(i) = temp(1);
      end
      Table(:,1) = start; 
      % 构建解空间
      points_index = 1:n;
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个点路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % 已访问的点集合(禁忌表)
             allow_index = ~ismember(points_index,tabu);    % 参加说明1（程序底部）
             allow = points_index(allow_index);  % 待访问的点集合
             P = allow;
             % 计算点之间转移概率
             for k = 1:length(allow)
                 P(k) = Tau(tabu(end),allow(k))^alpha * Heu_F(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % 轮盘赌法选择下一个访问的点,将整个P里的元素总和看作1，每个信息素元素表示一个扇形
            Pc = cumsum(P);   
            target_index = find(Pc >= rand);    %随机一个0到1之间的变量，找出信息素大于它的元素
            target = allow(target_index(1));    %轮盘赌法能保证每次循环都能随机命中概率较大的目标
            Table(i,j) = target;
         end
      end
      % 计算各个蚂蚁的路径距离
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % 计算最短路径距离及平均距离
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
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个点计算
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-vol) * Tau + Delta_Tau;
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
end

%% 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
Time_Cost=etime(clock,t0);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);
disp(['收敛迭代次数:' num2str(Limit_iter)]);
disp(['程序执行时间:' num2str(Time_Cost) '秒']);

%% 绘图
figure(1)
plot([points(Shortest_Route,1);points(Shortest_Route(1),1)],...  %三点省略符为Matlab续行符
     [points(Shortest_Route,2);points(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(points,1)
    text(points(i,1),points(i,2),['   ' num2str(i)]);
end
text(points(Shortest_Route(1),1),points(Shortest_Route(1),2),'       起点');
text(points(Shortest_Route(end),1),points(Shortest_Route(end),2),'       终点');
xlabel('点位置横坐标')
ylabel('点位置纵坐标')
title(['ACA最优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b')
legend('最短距离')
xlabel('迭代次数')
ylabel('距离')
title('算法收敛轨迹')
