
clear all,clc
ttt=clock;
	a = 0.99;	% �¶�˥�������Ĳ���
	t0 = 97; tf = 3; t = t0;
	Markov_length = 10000;	% Markov������
	points = xlsread('points_data.xlsx', 'B2:C53');
	numofpoi = size(points,1); 	% �����Ŀ
	% ͨ���������ķ�������������
	pointsidis = zeros(numofpoi, numofpoi);
	coor_x_tmp1 = points(:,1) * ones(1,numofpoi);
	coor_x_tmp2 = coor_x_tmp1';
	coor_y_tmp1 = points(:,2) * ones(1,numofpoi);
	coor_y_tmp2 = coor_y_tmp1';
	pointsidis = sqrt((coor_x_tmp1-coor_x_tmp2).^2 + ...
					(coor_y_tmp1-coor_y_tmp2).^2);
               xlswrite('dist',pointsidis)

	sol_new = 1:numofpoi;         % ������ʼ��
% sol_new��ÿ�β������½⣻sol_current�ǵ�ǰ�⣻sol_best����ȴ�е���ý⣻
	E_current = inf;E_best = inf; 		% E_current�ǵ�ǰ���Ӧ�Ļ�·���룻
% E_new���½�Ļ�·���룻
% E_best�����Ž��
	sol_current = sol_new; sol_best = sol_new;          
	p = 1;

	while t>=tf
		for r=1:Markov_length		% Markov������
			% ��������Ŷ�
			if (rand < 0.5)	% ��������ǽ�������������������
				% ������
				ind1 = 0; ind2 = 0;
				while (ind1 == ind2)
					ind1 = ceil(rand.*numofpoi);
					ind2 = ceil(rand.*numofpoi);
				end
				tmp1 = sol_new(ind1);
				sol_new(ind1) = sol_new(ind2);
				sol_new(ind2) = tmp1;
			else
				% ������
				ind1 = 0; ind2 = 0; ind3 = 0;
				while (ind1 == ind2) || (ind1 == ind3) ...
					|| (ind2 == ind3) || (abs(ind1-ind2) == 1)
					ind1 = ceil(rand.*numofpoi);
					ind2 = ceil(rand.*numofpoi);
					ind3 = ceil(rand.*numofpoi);
				end
				tmp1 = ind1;tmp2 = ind2;tmp3 = ind3;
				% ȷ��ind1 < ind2 < ind3
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

			%����Ƿ�����Լ��
			
			% ����Ŀ�꺯��ֵ�������ܣ�
			E_new = 0;
			for i = 1 : (numofpoi-1)
				E_new = E_new + ...
					pointsidis(sol_new(i),sol_new(i+1));
			end
			% �����ϴ����һ���㵽��һ����ľ���
			E_new = E_new + ...
				pointsidis(sol_new(numofpoi),sol_new(1));
			
			if E_new < E_current
				E_current = E_new;
				sol_current = sol_new;
				if E_new < E_best
% ����ȴ��������õĽⱣ������
					E_best = E_new;
					sol_best = sol_new;
				end
			else
				% ���½��Ŀ�꺯��ֵС�ڵ�ǰ��ģ�
				% �����һ�����ʽ����½�
				if rand < exp(-(E_new-E_current)./t)
					E_current = E_new;
					sol_current = sol_new;
				else	
					sol_new = sol_current;
				end
			end
		end
		t=t.*a;		% ���Ʋ���t���¶ȣ�����Ϊԭ����a��
    end
    
    Time_Cost=etime(clock,ttt);
    disp(['����ִ��ʱ��:' num2str(Time_Cost) '��']); 
    disp('���Ž�Ϊ��')
	disp(sol_best)
	disp('��̾��룺')
	disp(E_best)
    plot(points(:,1),points(:,2),'bo')
    axis([0 2000 0 1200]);
    hold on
    plot(points(sol_best,1),points(sol_best,2))
