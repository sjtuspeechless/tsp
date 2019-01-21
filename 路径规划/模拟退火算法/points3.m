
clear all,clc
ttt=clock;
	a = 0.99;	% 温度衰减函数的参数
	t0 = 97; tf = 3; t = t0;
	Markov_length = 10000;	% Markov链长度
	points = xlsread('points_data.xlsx', 'B2:C53');
	numofpoi = size(points,1); 	% 点的数目
	% 通过向量化的方法计算距离矩阵
	pointsidis = zeros(numofpoi, numofpoi);
	coor_x_tmp1 = points(:,1) * ones(1,numofpoi);
	coor_x_tmp2 = coor_x_tmp1';
	coor_y_tmp1 = points(:,2) * ones(1,numofpoi);
	coor_y_tmp2 = coor_y_tmp1';
	pointsidis = sqrt((coor_x_tmp1-coor_x_tmp2).^2 + ...
					(coor_y_tmp1-coor_y_tmp2).^2);
               xlswrite('dist',pointsidis)

	sol_new = 1:numofpoi;         % 产生初始解
% sol_new是每次产生的新解；sol_current是当前解；sol_best是冷却中的最好解；
	E_current = inf;E_best = inf; 		% E_current是当前解对应的回路距离；
% E_new是新解的回路距离；
% E_best是最优解的
	sol_current = sol_new; sol_best = sol_new;          
	p = 1;

	while t>=tf
		for r=1:Markov_length		% Markov链长度
			% 产生随机扰动
			if (rand < 0.5)	% 随机决定是进行两交换还是三交换
				% 两交换
				ind1 = 0; ind2 = 0;
				while (ind1 == ind2)
					ind1 = ceil(rand.*numofpoi);
					ind2 = ceil(rand.*numofpoi);
				end
				tmp1 = sol_new(ind1);
				sol_new(ind1) = sol_new(ind2);
				sol_new(ind2) = tmp1;
			else
				% 三交换
				ind1 = 0; ind2 = 0; ind3 = 0;
				while (ind1 == ind2) || (ind1 == ind3) ...
					|| (ind2 == ind3) || (abs(ind1-ind2) == 1)
					ind1 = ceil(rand.*numofpoi);
					ind2 = ceil(rand.*numofpoi);
					ind3 = ceil(rand.*numofpoi);
				end
				tmp1 = ind1;tmp2 = ind2;tmp3 = ind3;
				% 确保ind1 < ind2 < ind3
				if (ind1 < ind2) && (ind2 < ind3)
					;
				elseif (ind1 < ind3) && (ind3 < ind2)
					ind2 = tmp3;ind3 = tmp2;
				elseif (ind2 < ind1) && (ind1 < ind3)
					ind1 = tmp2;ind2 = tmp1;
				elseif (ind2 < ind3) && (ind3 < ind1) 
					ind1 = tmp2;ind2 = tmp3; ind3 = tmp1;
				elseif (ind3 < ind1) && (ind1 < ind2)
					ind1 = tmp3;ind2 = tmp1; ind3 = tmp2;
				elseif (ind3 < ind2) && (ind2 < ind1)
					ind1 = tmp3;ind2 = tmp2; ind3 = tmp1;
				end
				
				tmplist1 = sol_new((ind1+1):(ind2-1));
				sol_new((ind1+1):(ind1+ind3-ind2+1)) = ...
					sol_new((ind2):(ind3));
				sol_new((ind1+ind3-ind2+2):ind3) = ...
					tmplist1;
			end

			%检查是否满足约束
			
			% 计算目标函数值（即内能）
			E_new = 0;
			for i = 1 : (numofpoi-1)
				E_new = E_new + ...
					pointsidis(sol_new(i),sol_new(i+1));
			end
			% 再算上从最后一个点到第一个点的距离
			E_new = E_new + ...
				pointsidis(sol_new(numofpoi),sol_new(1));
			
			if E_new < E_current
				E_current = E_new;
				sol_current = sol_new;
				if E_new < E_best
% 把冷却过程中最好的解保存下来
					E_best = E_new;
					sol_best = sol_new;
				end
			else
				% 若新解的目标函数值小于当前解的，
				% 则仅以一定概率接受新解
				if rand < exp(-(E_new-E_current)./t)
					E_current = E_new;
					sol_current = sol_new;
				else	
					sol_new = sol_current;
				end
			end
		end
		t=t.*a;		% 控制参数t（温度）减少为原来的a倍
    end
    
    Time_Cost=etime(clock,ttt);
    disp(['程序执行时间:' num2str(Time_Cost) '秒']); 
    disp('最优解为：')
	disp(sol_best)
	disp('最短距离：')
	disp(E_best)
    plot(points(:,1),points(:,2),'bo')
    axis([0 2000 0 1200]);
    hold on
    plot(points(sol_best,1),points(sol_best,2))
