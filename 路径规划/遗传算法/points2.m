 
 %% 加载问题的数据
 ttt=clock;
points =xlsread('points_data.xlsx', 'B2:C53');
numofpo=size(points,1);
plot(points(:,1),points(:,2),'bo');
xlabel('点的横坐标x'); ylabel('点的纵坐标y'); 
grid on

%% 计算点之间距离

pointdis = zeros(numofpo);
for count1=1:numofpo
    for count2=1:count1
        x1 = points(count1,1);
        y1 = points(count1,2);
        x2 = points(count2,1);
        y2 = points(count2,2);
        pointdis(count1,count2)=sqrt((x1-x2)^2+(y1-y2)^2);
        pointdis(count2,count1)=pointdis(count1,count2);
    end
end


%% 定义目标函数
FitnessFcn = @(x) points_fitness(x,pointdis);   %适应度

my_plot = @(options,state,flag) points_plot(options, ...
    state,flag,points);

%% 设置优化属性并执行遗传算法求解

options = optimoptions(@ga, 'PopulationType', 'custom','InitialPopulationRange', ...
                            [1;numofpo]);       %创建优化选项

options = optimoptions(options,'CreationFcn',@create_permutations, ...  %初始化填充函数
                        'CrossoverFcn',@crossover_permutation, ...      %创建交叉子集的函数
                        'MutationFcn',@mutate_permutation, ...          %产生变异后代的函数
                        'PlotFcn', my_plot, ...                         %用于绘制由算法计算的数据的函数
                        'MaxGenerations',1000,'PopulationSize',100, ...   %设置人口大小
                        'MaxStallGenerations',200,'UseVectorized',true);

numberOfVariables = numofpo;
[x,fval,reason,output] = ...
    ga(FitnessFcn,numberOfVariables,[],[],[],[],[],[],[],options)
Time_Cost=etime(clock,ttt);
    disp(['程序执行时间:' num2str(Time_Cost) '秒']); 
