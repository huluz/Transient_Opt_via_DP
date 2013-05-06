%采用隐式算法进行非稳态模拟

clc;clear;

%气体参数
mol = 17.44;            		%气体分子质量
R = 8314;               		%通用气体常数
PreC = 46.59E5;         		%临界压力
TempC = 198.19;         		%临界温度
Temp = 15+273.15;       	%输气温度
Temp_sta = 20+273.15;   	%标况温度
Pre_sta = 101325;       		%标况压力
Den_air = 1.209;        		%空气密度
Rel_Den = 0.6;          		%天然气相对密度
Den_sta = Rel_Den*Den_air;  	%天然气标况密度
alpha = Den_sta*Temp_sta/Pre_sta/Temp;  	%通用气体方程参数
beta = 0.257/PreC - 0.533*TempC/PreC/Temp; 	%AGA气体状态方程参数

%管段及模拟参数
Len = 100E3;            		%管段长度
lamda = 0.008;          		%摩阻系数
Din = 0.6096;           		%管段内径
C0 = 0.03848;           		%稳态模拟公式系数
Time = 24*3600;		%模拟时长
Time_Sec = 3600;		%单个时间段时长
Secs = Time/Time_Sec;	%时间段数目
dt = 60 * 15;              		%时步
TimeSteps_Total = Time / dt; 	%总时步数
TimeSteps_Per_Sec = Time_Sec / dt;		%每个时间段的步长数
dx = 10E3;               		%空间步长
SpaceSteps = Len / dx;		%空间分段数
Area = 0.25*pi*Din^2;		%管段横截面积
Volumn_Sec = Area*dx;
gama = 1;			%阻尼系数

%边界条件
Storage = zeros(TimeSteps_Total,1);		%管存量
Qbasic = 65;           		 		%流量基数
Qs_Opt = [25; 5; 5; 5; 45; 45; 105; 65; 45; 85; 85; 125; 125; 5; 85; 85; 125; 145; 145; 65; 45; 65; 5; 25];	%起点流量
Qs = zeros(TimeSteps_Total,1);
for ii = 1:Secs 					%根据时间点上的值设定整个时间段的流量
	%Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
	for j = 1:TimeSteps_Per_Sec
		if ii == 1
			Qs(TimeSteps_Per_Sec*(ii-1)+j) = (Qs_Opt(ii)-Qs_Opt(Secs))*j/TimeSteps_Per_Sec + Qs_Opt(Secs);
		else
			Qs(TimeSteps_Per_Sec*(ii-1)+j) = (Qs_Opt(ii)-Qs_Opt(ii-1))*j/TimeSteps_Per_Sec + Qs_Opt(ii-1);
		end
	end
end
Mss = (Den_sta/Area)*Qs;			%起点质量流量密度
Ps = zeros(TimeSteps_Total,1);		%起点压力
Pe = zeros(TimeSteps_Total,1);		%终点压力
Ff = [0.2; 0.15; 0.1; 0.25; 0.35; 0.58; 1.2; 1.3; 1; 0.97; 0.85; 1.65; 2; 1; 0.8; 0.65; ...
    1.15; 1.9; 2.8; 2.2; 1.2; 0.85; 0.5; 0.35];    		%小时流量不均匀系数
Qe = zeros(TimeSteps_Total,1);         		%终点流量
for i = 1:Secs 					%根据时间点上的值设定整个时间段的流量
	%Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
	for j = 1:TimeSteps_Per_Sec
		if i == 1
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(Secs))*j/TimeSteps_Per_Sec + Ff(Secs));
		else
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(i-1))*j/TimeSteps_Per_Sec + Ff(i-1));
		end
	end
end
Mse = Den_sta * Qe/Area;			%终点质量流量密度

%稳态模拟
tl = Len;               		%管段长度
Pls = 5.5e6;            		%起点压力
i = 1;
Pressure(i) = Pls;		%沿线压力记录
while tl>0              		%稳态模拟
    z = 1 + beta*Pls;   		%压缩因子
    Ple = Pls^2 - lamda*z*Rel_Den*Temp*dx*Qe(TimeSteps_Total)^2/C0^2/Din^5;
    Ple = Ple^0.5;
    i = i+1;
    Pressure(i) = Ple;
    tl = tl - dx;
    Pls = Ple;
end
MassFlux = (Den_sta*Qe(TimeSteps_Total)/Area)*ones(SpaceSteps+1,1);	%构造初始条件-质量流量密度
Pressure = Pressure';		%压力
Storage_Total = 0;		%计算管存量
for ii = 1:SpaceSteps
	Pre_Aver = 2*(Pressure(ii) + Pressure(ii + 1)^2/(Pressure(ii) + Pressure(ii + 1)))/3;
	z = 1 + beta*Pre_Aver;
	Den_Rel = Pre_Aver/z/R/Temp;
	Storage_Sec = Volumn_Sec*Den_Rel;
	Storage_Total = Storage_Total + Storage_Sec;
end
disp(['Initial Total Storage: ' sprintf('%d',Storage_Total)]);
%figure(1);
%plot(Pressure);
%title('Pressure');
%figure(2);
%plot(MassFlux);
%title('MassFlux');

%瞬态模拟
%构造管流控制方程组
creat_transfun_re01(SpaceSteps);
tic;
for i = 1:TimeSteps_Total
	tf = @(x)transfun_re01(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Mss(i),Mse(i));	%构造方程
	x0 = zeros(2*SpaceSteps,1);	%准备初值
	x0(1) = Pressure(1);
	for j = 2:SpaceSteps
		x0(2*j-2) = Pressure(j);
		x0(2*j-1) = MassFlux(j);
	end
	x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
	%options = optimset('MaxFunEvals',1E4,'MaxIter',1E4,'Display','iter-detailed');
	options = optimset('Display','off');
	results = fsolve(tf,x0,options);		%计算
	results = (1-gama)*x0+gama*results;	%引入阻尼系数，防止计算结果震荡
	Pressure(1) = results(1);	%归档计算结果
	for j = 2: SpaceSteps
		Pressure(j) = results(2*j-2);
		MassFlux(j) = results(2*j-1);
	end
	Pressure(SpaceSteps+1) = results(2*SpaceSteps);
	Ps(i) = results(1);		%管段起点质量流量
	Pe(i) = results(2*SpaceSteps);	%管段终点压力
	if i < TimeSteps_Total
		MassFlux(1) = Mss(i+1);	%引入边界条件
		MassFlux(SpaceSteps+1) = Mse(i+1);
	end
	Storage_Total = 0;		%计算管存量
	for ii = 1:SpaceSteps
		Pre_Aver = 2*(Pressure(ii) + Pressure(ii + 1)^2/(Pressure(ii) + Pressure(ii + 1)))/3;
		z = 1 + beta*Pre_Aver;
		Den_Rel = Pre_Aver/z/R/Temp;
		Storage_Sec = Volumn_Sec*Den_Rel;
		Storage_Total = Storage_Total + Storage_Sec;
	end
	Storage(i) = Storage_Total;
	if mod(i,10) == 0
		%计算过程图形化
		figure(3);			%终点压力变化
		plot(Pe(1:i));
		title('Pressure at the Outlet');
		figure(4);			%起点压力变化
		plot(Ps(1:i));
		title('Pressure at the Inlet');
		figure(5);			%管存量
		plot(Storage(1:i));
		title('Storage of the Pipeleg');
	end
end
toc;

%计算结果可视化
figure(3);			%终点压力变化
plot(Pe);
title('Pressure at the Outlet');
Qs = (Area/Den_sta)*Mss;	%起点压力变化
figure(4);
plot(Ps);
title('Pressure at the Inlet');
figure(5);			%管存量
plot(Storage);
title('Storage of the Pipeleg');