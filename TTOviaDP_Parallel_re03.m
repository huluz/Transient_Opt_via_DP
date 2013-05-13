function TTOviaDP_Parallel_re03(Ff, Time, dt, Qs_Min, Qs_Max, dq, Qbasic, FN)
	% 采用DP算法求解输气管道末段非稳态优化问题
	%Ff - 小时不均匀用气系数
	%Time - 优化时间段，s
	%dt - 时步，s
	%Qs_Min, Qs_Max - 管段进口流量范围, Nm^3/s
	%dq - 流量离散步长, Nm^3/s
	%Qbasic - 基础流量, Nm^3/s
	%FN - 计算结果存储文件名

%气体相关参数
mol = 17.44;			%相对分子质量
R = 8314;			%通用气体常数
PreC = 46.59E5;		%气体视临界压力
TempC = 198.19;		%气体视临界温度
Temp = 15 + 273.15;		%输气温度
Den_air = 1.209;		%标况下空气密度，kg/m^3
Temp_sta = 20 + 273.15;	%标况温度
Pre_sta = 101325;		%标况压力
Rel_Den = 0.6;			%相对密度
Den_sta = Rel_Den * Den_air;	%标况气体密度
alpha = Den_sta*Temp_sta/Pre_sta/Temp;	%通用气体状态方程参数
beta = 0.257/PreC - 0.533*TempC/PreC/Temp;	%AGA气体状态方程参数

%管段参数
Len = 100E3;			%管道长度，m
Pin = 3.5E6;			%压缩机入口压力
lamda = 0.008;			%摩阻系数
Din = 0.6096;			%管内径，m
Area = 0.25*pi*Din^2;		%管段横截面积
Pe_min = 4.5E6;		%管段出口允许最低压力
%Ps_max = 7.3E6;		%管段进口允许最高压力

%计算参数
Time_Per_Sec = 3600;		%流量边界条件设定时步，s
Time_Secs = Time / Time_Per_Sec;	%时间段数
TimeSteps_Total = Time / dt;	%总时步数
TimeSteps_Per_Sec = Time_Per_Sec / dt;	%每个时间段包含的时步数
%dp = 0.2E6;			%状态变量空间离散步长
dx = 10E3;			%空间步长，m
SpaceSteps = Len / dx;		%空间分段数
creat_transfun_re01(SpaceSteps);		%创建状态转移方程

%模拟参数
C0 = 0.03848;			%稳态模拟公式参数
%a = [4.4889;0.071;-0.009];	%压缩机压头-流量曲线系数
%b = [9.0719;2.2055;-0.0161];	%压缩机效率-流量曲线系数
Pressure_ini = zeros(SpaceSteps+1,1);		%沿线压力分布
MassFlux_ini = zeros(SpaceSteps+1,1);	%沿线质量流量分布

%优化计算参数

%边界条件
Qe = zeros(TimeSteps_Total + 1,1);		%终点流量
Qe(1) = Ff(1)*Qbasic;
for i = 1:Time_Secs			%根据时间点上的值设定整个时间段的流量
	for j = 1:TimeSteps_Per_Sec
		Qe(TimeSteps_Per_Sec*(i-1)+j+1) = Qbasic*((Ff(i+1)-Ff(i))*j/TimeSteps_Per_Sec + Ff(i));
	end
end
Mse = Den_sta * Qe/Area;	%终点质量流量密度
%figure;plot(Ff);figure;plot(Qe);

%初始条件
tl = Len;			%管段长度
Ple = Pe_min;			%起点压力
Pressure_ini(SpaceSteps+1) = Ple;%沿线压力记录
i = SpaceSteps;
while tl>0 			%稳态模拟
	z = 1 + beta*Ple;	%压缩因子
	Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(1)^2/C0^2/Din^5;
	Pls = Pls^0.5;
	Pressure_ini(i) = Pls;
	i = i - 1;
	tl = tl - dx;
	Ple = Pls;
end
MassFlux_ini = (Den_sta*Qe(1)/Area)*ones(SpaceSteps+1,1);	%构造初始条件-质量流量密度
Pressure_ini = Pressure_ini';		%压力
%figure(1);
%plot(Pressure_ini);
%title('Pressure_ini');
%figure(2);
%plot(MassFlux_ini);
%title('MassFlux_ini');

%DP算法求解优化问题

%构建DP算法初始状态
Results_Now = [Pressure_ini,MassFlux_ini',zeros(1,Time_Secs),0];	%设定初始状态							%管存量
Volumn_Sec = Area*dx;
Col_Num_Per_Rec = 2*(SpaceSteps + 1) + Time_Secs + 1;		%每个记录所需的行数
Ms_Rec_Start_Num = SpaceSteps + 2;					%质量流量记录起始位置
Desicision_Rec_Start_Num = 2*(SpaceSteps + 1) + 1;			%决策序列记录起始问题

%离散可行域
Qs_avai_Num = (Qs_Max - Qs_Min) / dq + 1;		%离散点数量
Qs_avai = zeros(Qs_avai_Num, 1);			%确定离散点数值
for ii = 1:Qs_avai_Num
	Qs_avai(ii) = Qs_Min + (ii - 1)*dq;
end
Qs_avai = (Den_sta/Area)*Qs_avai;			%转换单位

%顺序递推过程
for i = 2:Time_Secs+1
	startmatlabpool;
	disp('=========================================');
	disp(['Time: ' sprintf('%d',i - 1)]);
	tic;				%计时

	%转储上一阶段的计算结果，并准备继续计算所需的内存空间
	Results_Pre = Results_Now;			%转储上一阶段的计算结果
	States_Pre_Num = size(Results_Pre,1);	%本阶段/上一时步状态数
	States_Now_Num_Temp = States_Pre_Num*Qs_avai_Num;			%下一阶段/本时步状态数初步估计
	Results_Now_Temp = zeros(States_Now_Num_Temp,Col_Num_Per_Rec);	%准备计算所需内存空间
	Results_Now_Temp_Pressure = zeros(States_Now_Num_Temp, Ms_Rec_Start_Num - 1);				%压力记录
	Results_Now_Temp_MassFlux = zeros(States_Now_Num_Temp, Desicision_Rec_Start_Num - Ms_Rec_Start_Num);	%质量流量记录
	Results_Now_Temp_Decision = zeros(States_Now_Num_Temp, Col_Num_Per_Rec - Desicision_Rec_Start_Num);	%决策过程记录
	Results_Now_Temp_Decision_TT = zeros(States_Now_Num_Temp,1);
	Results_Now_Temp_Com_Consum = zeros(States_Now_Num_Temp, 1);						%压缩机能耗记录
	Storage = zeros(States_Now_Num_Temp,1);	%管存量

	%产生状态与决策变量组合，并进行模拟
	parfor Rec_Num = 1:States_Now_Num_Temp 		%按照估算状态数进行循环
		%解析循环变量，得到状态与决策变量索引
		n = mod(Rec_Num, Qs_avai_Num) ;		%决策变量索引
		if n == 0 
			n = Qs_avai_Num;
			m = Rec_Num/Qs_avai_Num;
		else
			m = (Rec_Num - n)/Qs_avai_Num + 1;	%状态索引
		end

		%复制决策过程
		Results_Now_Temp_Decision(Rec_Num, :) = Results_Pre(m, Desicision_Rec_Start_Num:Col_Num_Per_Rec-1);
		Results_Now_Temp_Decision_TT(Rec_Num) = Qs_avai(n);				%设定决策变量值
		Pressure = Results_Pre(m,1:Ms_Rec_Start_Num - 1);					%设定压力初始值，准备模拟
		MassFlux = Results_Pre(m,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1);	%设定质量流量初始值

		%模拟过程
		Qs_Pre = Results_Pre(m,Ms_Rec_Start_Num);	%取出上一时步管段入口处流量
		Com_Consum = 0;					%压缩机能耗
		%单个时间段内可能存在多个时步，需要进行多次模拟
		for l = 1:TimeSteps_Per_Sec
			Qs_Simu = (Qs_avai(n) - Qs_Pre)*l/TimeSteps_Per_Sec + Qs_Pre;	%通过插值方法构建边界条件
			tf = @(x)transfun_re01(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Qs_Simu,Mse((i-2)*TimeSteps_Per_Sec+l+1));	%构造方程
			x0 = zeros(2*SpaceSteps,1);	%准备初值
			x0(1) = Pressure(1);
			for j = 2:SpaceSteps
				x0(2*j-2) = Pressure(j);
				x0(2*j-1) = MassFlux(j);
			end
			x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
			options = optimset('Display','off');
			results = fsolve(tf,x0,options);	%计算

			Pressure(1) = results(1);	%归档计算结果
			for j = 2:SpaceSteps
				Pressure(j) = results(2*j-2);
				MassFlux(j) = results(2*j-1);
			end
			Pressure(SpaceSteps+1) = results(2*SpaceSteps);
			MassFlux(1) = Qs_Simu;
			MassFlux(SpaceSteps+1) = Mse((i-2)*TimeSteps_Per_Sec+l+1);

			%计算压缩机能耗
			if Pressure(SpaceSteps+1) < Pe_min - 0.1e6 && i == 2 					%初始时放宽一些约束条件
				Results_Now_Temp_Com_Consum(Rec_Num) = -1;
				break;
			elseif  Pressure(SpaceSteps + 1) < Pe_min 						%后续过程收紧约束条件
				Results_Now_Temp_Com_Consum(Rec_Num) = -1;
				break;
			else
				Quan_Temp = (0.328*Area/Den_sta)*MassFlux(1);
				Sec_Com_Consum = dt*Quan_Temp*(2.682*(Pressure(1)/Pin)^0.217 - 2.658);	%该时步压缩机功率
				Com_Consum = Com_Consum + Sec_Com_Consum;
			end
		end
		%归档计算结果
		if Results_Now_Temp_Com_Consum(Rec_Num) ~= -1
			Results_Now_Temp_Com_Consum(Rec_Num) = Results_Pre(m,Col_Num_Per_Rec) + Com_Consum;
			for ii = 1:SpaceSteps
				Pre_Aver = 2*(Pressure(ii) + Pressure(ii + 1)^2/(Pressure(ii) + Pressure(ii + 1)))/3;
				z = 1 + beta*Pre_Aver;
				Den_Rel = Pre_Aver/z/R/Temp;
				Storage_Sec = Volumn_Sec*Den_Rel;
				Storage(Rec_Num) = Storage(Rec_Num) + Storage_Sec;
			end
		end
		Results_Now_Temp_Pressure(Rec_Num,:) = Pressure;
		Results_Now_Temp_MassFlux(Rec_Num,:) = MassFlux;
	end

	%整理计算结果
	for Rec_Num = 1:States_Now_Num_Temp
		Results_Now_Temp(Rec_Num,1:Ms_Rec_Start_Num - 1) = Results_Now_Temp_Pressure(Rec_Num,:);				%压力记录
		Results_Now_Temp(Rec_Num,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1) = Results_Now_Temp_MassFlux(Rec_Num,:);	%质量流量记录
		Results_Now_Temp_Decision(Rec_Num,i-1) = Results_Now_Temp_Decision_TT(Rec_Num);						%决策过程记录
		Results_Now_Temp(Rec_Num,Desicision_Rec_Start_Num:Col_Num_Per_Rec-1) = Results_Now_Temp_Decision(Rec_Num,:);
		Results_Now_Temp(Rec_Num,Col_Num_Per_Rec) = Results_Now_Temp_Com_Consum(Rec_Num);					%能耗记录
	end

	%根据管存收缩状态空间
	Red_States = 0;	%冗余状态数
	for kk = 1:States_Now_Num_Temp - 1
		if Results_Now_Temp(kk,Col_Num_Per_Rec) ~= -1
			for mm = kk+1 : States_Now_Num_Temp
				if Results_Now_Temp(mm,Col_Num_Per_Rec) ~= -1 && Storage(kk) > Storage(mm) && Results_Now_Temp(kk,Col_Num_Per_Rec) < Results_Now_Temp(mm,Col_Num_Per_Rec)
					Results_Now_Temp(mm, Col_Num_Per_Rec) = -1;
					Red_States = Red_States + 1;
				end
			end
		end
	end

	%剔除不可行方案
	Bad_Recs = 0;		%统计可行方案数目
	for l = 1:States_Now_Num_Temp
		if Results_Now_Temp(l,Col_Num_Per_Rec) == -1
			Bad_Recs = Bad_Recs + 1;
		end
	end
	Good_Recs = States_Now_Num_Temp - Bad_Recs;
	if Good_Recs <= 0
		sprintf('%s\n','No available results!');
		return;
	end
	l = 1;m=1;		%归档可行方案
	Results_Now = zeros(Good_Recs,Col_Num_Per_Rec);
	while l <= Good_Recs
		if Results_Now_Temp(m,Col_Num_Per_Rec) ~= -1
			Results_Now(l,:) = Results_Now_Temp(m,:);
			l = l + 1;
		end
		m = m + 1;
	end

	%计算过程可视化
	disp(['Total Combinations: ' sprintf('%d',States_Now_Num_Temp)]);
	disp(['Good Recodes: ' sprintf('%d',Good_Recs)]);
	disp(['Bad Recodes: ' sprintf('%d',Bad_Recs)]);
	disp(['Redundancy Recodes: ' sprintf('%d', Red_States)]);
	toc;			%单个时间段计算时间
	disp('=========================================');
	disp('');
end

%筛选最优方案
%Rec_Num = size(Results_Now,1);
%Min_Obj = Results_Now(1,Col_Num_Per_Rec);
%Opt_Rec_Num = 1;
%for i = 1:Rec_Num
%	if Min_Obj > Results_Now(i,Col_Num_Per_Rec)
%		Min_Obj = Results_Now(i,Col_Num_Per_Rec);
%		Opt_Rec_Num = i;
%	end
%end
[Min_Obj, Opt_Rec_Num] = min(Results_Now(:,Col_Num_Per_Rec));
Opt_Des = Area/Den_sta*Results_Now(Opt_Rec_Num,Desicision_Rec_Start_Num:Col_Num_Per_Rec-1);

%图形化计算结果
%figure(3);
%plot(Opt_Des);
%title('Optimum Control Strategy');
%legend('Volumn quantity at the inlet');

save(FN);