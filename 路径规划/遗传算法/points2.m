 
 %% �������������
 ttt=clock;
points =xlsread('points_data.xlsx', 'B2:C53');
numofpo=size(points,1);
plot(points(:,1),points(:,2),'bo');
xlabel('��ĺ�����x'); ylabel('���������y'); 
grid on

%% �����֮�����

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


%% ����Ŀ�꺯��
FitnessFcn = @(x) points_fitness(x,pointdis);   %��Ӧ��

my_plot = @(options,state,flag) points_plot(options, ...
    state,flag,points);

%% �����Ż����Բ�ִ���Ŵ��㷨���

options = optimoptions(@ga, 'PopulationType', 'custom','InitialPopulationRange', ...
                            [1;numofpo]);       %�����Ż�ѡ��

options = optimoptions(options,'CreationFcn',@create_permutations, ...  %��ʼ����亯��
                        'CrossoverFcn',@crossover_permutation, ...      %���������Ӽ��ĺ���
                        'MutationFcn',@mutate_permutation, ...          %�����������ĺ���
                        'PlotFcn', my_plot, ...                         %���ڻ������㷨��������ݵĺ���
                        'MaxGenerations',1000,'PopulationSize',100, ...   %�����˿ڴ�С
                        'MaxStallGenerations',200,'UseVectorized',true);

numberOfVariables = numofpo;
[x,fval,reason,output] = ...
    ga(FitnessFcn,numberOfVariables,[],[],[],[],[],[],[],options)
Time_Cost=etime(clock,ttt);
    disp(['����ִ��ʱ��:' num2str(Time_Cost) '��']); 
