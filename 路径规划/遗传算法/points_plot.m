function state = points_plot(options,state,flag,points)
if strcmpi(flag,'init')
  points =xlsread('points_data.xlsx', 'B2:C53');
end
[~,i] = min(state.Score);
genotype = state.Population{i};

plot(points(:,1),points(:,2),'bo');

hold on
plot(points(genotype,1),points(genotype,2));
xlabel('点的横坐标x'); ylabel('点的纵坐标y'); 
grid on
hold off
