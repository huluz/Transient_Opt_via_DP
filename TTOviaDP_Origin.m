% 采用DP算法求解输气管道末段非稳态优化问题

% 代码初始化
clc;clear;

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
Ps_max = 7.3E6;		%管段进口允许最高压力
Qbasic = 65;			%基础流量，m^3/s
Ff = [0.2; 0.15; 0.1; 0.25; 0.35; 0.58; 1.2; 1.3; 1; 0.97; 0.85; 1.65; 2; 1; 0.8; 0.65; ...
	1.15; 1.9; 2.8; 2.2; 1.2; 0.85; 0.5; 0.35];	%小时不均匀用气系数

%计算参数
Time = 24*3600;		%优化时间段，s
Time_Per_Sec = 3600;		%流量边界条件设定时步，s
dt = 1800;			%时步，s
Time_Secs = Time / Time_Per_Sec;	%时间段数
TimeSteps_Total = Time / dt;	%总时步数
TimeSteps_Per_Sec = Time_Per_Sec / dt;	%每个时间段包含的时步数
dp = 0.2E6;			%状态变量空间离散步长
dx = 10E3;			%空间步长，m
SpaceSteps = Len / dx;		%空间分段数
creat_transfun(SpaceSteps);	%创建状态转移方程

%模拟参数
C0 = 0.03848;			%稳态模拟公式参数
a = [4.4889;0.071;-0.009];	%压缩机压头-流量曲线系数
b = [9.0719;2.2055;-0.0161];	%压缩机效率-流量曲线系数
Pressure_ini = zeros(SpaceSteps+1,1);		%沿线压力分布
MassFlux_ini = zeros(SpaceSteps+1,1);	%沿线质量流量分布

%优化计算参数

%边界条件
Qe = zeros(TimeSteps_Total,1);		%终点流量
for i = 1:Time_Secs			%根据时间点上的值设定整个时间段的流量
	%Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
	for j = 1:TimeSteps_Per_Sec
		if i == 1
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(Time_Secs))*j/TimeSteps_Per_Sec + Ff(Time_Secs));
		else
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(i-1))*j/TimeSteps_Per_Sec + Ff(i-1));
		end
	end
end
Mse = Den_sta * Qe/Area;	%终点质量流量密度

%DP算法初始化-设置初始状态
%稳态模拟
tl = Len;			%管段长度
Ple = Pe_min;			%起点压力
Pressure_ini(SpaceSteps+1) = Ple;%沿线压力记录
i = SpaceSteps;
while tl>0 			%稳态模拟
	z = 1 + beta*Ple;	%压缩因子
	Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(TimeSteps_Total)^2/C0^2/Din^5;
	Pls = Pls^0.5;
	Pressure_ini(i) = Pls;
	i = i - 1;
	tl = tl - dx;
	Ple = Pls;
end
MassFlux_ini = (Den_sta*Qe(TimeSteps_Total)/Area)*ones(SpaceSteps+1,1);	%构造初始条件-质量流量密度
Pressure_ini = Pressure_ini';		%压力
%figure(1);
%plot(Pressure_ini);
%title('Pressure_ini');
%figure(2);
%plot(MassFlux_ini);
%title('MassFlux_ini');

%DP算法求解优化问题

%构建DP算法初始状态
Results_Now = [Pressure_ini,MassFlux_ini',zeros(1,Time_Secs),0];	%设定初始状态
%Storage = 0;								%管存量
Volumn_Sec = Area*dx;
%for ii = 1:SpaceSteps
%	Pre_Aver = (Pressure_ini(ii) + Pressure_ini(ii+1))/2;
%	z = 1 + beta*Pre_Aver;
%	Den_Rel = Pre_Aver/z/R/Temp;
%	Storage_Sec = Volumn_Sec*Den_Rel;
%	Storage = Storage + Storage_Sec;
%end
Col_Num_Per_Rec = 2*(SpaceSteps + 1) + Time_Secs + 1;		%每个记录所需的行数
Ms_Rec_Start_Num = SpaceSteps + 2;					%质量流量记录起始位置
Desicision_Rec_Start_Num = 2*(SpaceSteps + 1) + 1;			%决策序列记录起始问题

%顺序递推过程
for i = 1:2
	disp('=========================================');
	disp(['Time: ' sprintf('%d',i)]);
	tic;				%计时
	%确定决策变量范围
	tl = Len;			%管段长度
	Ple = Pe_min;			%起点压力
	while tl>0 			%稳态模拟
		z = 1 + beta*Ple;	%压缩因子
		Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(i*TimeSteps_Per_Sec)^2/C0^2/Din^5;
		Pls = Pls^0.5;
		tl = tl - dx;
		Ple = Pls;
	end
	Ps_min = Pls;

	%离散决策变量可行域
	Ps_avai_Num = 0; 		%计算离散后得到的决策变量可取值数目
	while Ps_max - Ps_avai_Num*dp > Ps_min
		Ps_avai_Num = Ps_avai_Num + 1;
	end
	if Ps_avai_Num <= 0
		sprintf('%s\n','Desicision is Empty!');
		return;
	end
	%Ps_avai_Num
	Ps_avai = zeros(Ps_avai_Num,1);	%计算决策变量可取值
	for j = 1:Ps_avai_Num
		Ps_avai(j) = Ps_max - (j - 1)*dp;
	end

	%转储上一阶段的计算结果，并准备继续计算所需的内存空间
	Results_Pre = Results_Now;			%转储上一阶段的计算结果
	States_Pre_Num = size(Results_Pre,1);	%本阶段/上一时步状态数
	States_Now_Num_Temp = States_Pre_Num*Ps_avai_Num;			%下一阶段/本时步状态数初步估计
	Results_Now_Temp = zeros(States_Now_Num_Temp,Col_Num_Per_Rec);	%准备计算所需内存空间
	%Storage = zeros(States_Now_Num_Temp,1);	%管存量

	%产生状态与决策变量组合，并进行模拟
	%startmatlabpool
	for m = 1:States_Pre_Num 	%循环设定状态
		for n = 1:Ps_avai_Num 	%循环设定决策变量
			%为方便计算，准备一些参数
			Rec_Num = (m - 1)*Ps_avai_Num + n;		%当前记录号

			%复制决策过程
			Results_Now_Temp(Rec_Num, Desicision_Rec_Start_Num:Col_Num_Per_Rec - 1) = Results_Pre(m,Desicision_Rec_Start_Num:Col_Num_Per_Rec - 1);
			Results_Now_Temp(Rec_Num, Desicision_Rec_Start_Num + i - 1) = Ps_avai(n);	%设定决策变量值
			Pressure = Results_Pre(m,1:Ms_Rec_Start_Num - 1);					%设定压力初始值，准备模拟
			MassFlux = Results_Pre(m,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1);	%设定质量流量初始值

			%模拟过程
			Ps_Pre = Results_Pre(m,1);	%取出上一时步管段入口处压力
			Com_Consum = 0;		%压缩机能耗
			%单个时间段内可能存在多个时步，需要进行多次模拟
			for l = 1:TimeSteps_Per_Sec
				Ps_Simu = (Ps_avai(n) - Ps_Pre)*l/Time_Per_Sec + Ps_Pre;	%通过插值方法构建边界条件
				tf = @(x)transfun(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Ps_Simu,Mse((i-1)*TimeSteps_Per_Sec+l));	%构造方程
				x0 = zeros(2*SpaceSteps,1);	%准备初值
				x0(1) = MassFlux(1);
				for j = 2:SpaceSteps
					x0(2*j-2) = Pressure(j);
					x0(2*j-1) = MassFlux(j);
				end
				x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
				options = optimset('Display','off');
				results = fsolve(tf,x0,options);	%计算

				MassFlux(1) = results(1);	%归档计算结果
				for j = 2:SpaceSteps
					Pressure(j) = results(2*j-2);
					MassFlux(j) = results(2*j-1);
				end
				Pressure(SpaceSteps+1) = results(2*SpaceSteps);
				Pressure(1) = Ps_Simu;
				MassFlux(SpaceSteps+1) = Mse((i-1)*TimeSteps_Per_Sec+l);

				%计算压缩机能耗
				if Pressure(SpaceSteps+1) < Pe_min
					Results_Now_Temp(Rec_Num, Col_Num_Per_Rec) = -1;
					Com_Consum = 0;
				else
					Com_Consum = Com_Consum + dt*MassFlux(1)*Area*(Pressure(1)/Pin)^0.8;
				end
			end
			%归档计算结果
			if Results_Now_Temp(Rec_Num, Col_Num_Per_Rec) ~= -1
				Results_Now_Temp(Rec_Num, Col_Num_Per_Rec) = Results_Pre(m,Col_Num_Per_Rec) + Com_Consum;
%				for ii = 1:SpaceSteps
%					Pre_Aver = (Pressure(ii) + Pressure(ii+1))/2;
%					z = 1 + beta*Pre_Aver;
%					Den_Rel = Pre_Aver/z/R/Temp;
%					Storage_Sec = Volumn_Sec*Den_Rel;
%					Storage(Rec_Num) = Storage(Rec_Num) + Storage_Sec;
%				end
			end
			Results_Now_Temp(Rec_Num,1:Ms_Rec_Start_Num - 1) = Pressure;
			Results_Now_Temp(Rec_Num,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1) = MassFlux;
		end
	end

	%根据管存收缩状态空间
%	Red_States = 0;	%冗余状态数
%	for kk = 1:States_Now_Num_Temp - 1
%		if Results_Now_Temp(kk,Col_Num_Per_Rec) ~= -1
%			for mm = kk+1 : States_Now_Num_Temp
%				if Results_Now_Temp(mm,Col_Num_Per_Rec) ~= -1 && Storage(kk) > Storage(mm) && Results_Now_Temp(kk,Col_Num_Per_Rec) < Results_Now_Temp(mm,Col_Num_Per_Rec)
%					Results_Now_Temp(mm, Col_Num_Per_Rec) = -1;
%					Red_States = Red_States + 1;
%				end
%			end
%		end
%	end

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
	disp(['Bad Recodes; ' sprintf('%d',Bad_Recs)]);
	%disp(['Redundancy Recodes: ' sprintf('%d', Red_States)]);
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
Opt_Des = Results_Now(Opt_Rec_Num,Desicision_Rec_Start_Num:Col_Num_Per_Rec-1);

%图形化计算结果
%figure(3);
%plot(Opt_Des);
%title('Optimum Control Strategy');
%legend('Volumn quantity at the inlet');