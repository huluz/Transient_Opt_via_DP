function SimuIm_re03(fn)
	%fn - 存储有进行模拟计算边界条件的数据文件

%采用隐式算法进行非稳态模拟
load (fn, 'Time', 'Ff', 'Opt_Des', 'Qbasic','Pin');	%由计算结果载入模拟所需边界条件
Qs_Opt = [Ff(1)*Qbasic;Opt_Des'];

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
Time_Sec = 3600;		%单个时间段时长
Secs = Time/Time_Sec;	%时间段数目
dt = 60 * 10;              		%时步
TimeSteps_Total = Time / dt; 	%总时步数
TimeSteps_Per_Sec = Time_Sec / dt;		%每个时间段的步长数
dx = 10E3;               		%空间步长
SpaceSteps = Len / dx;		%空间分段数
Area = 0.25*pi*Din^2;		%管段横截面积
Volumn_Sec = Area*dx;
gama = 1;			%阻尼系数

%边界条件
Storage = zeros(TimeSteps_Total,1);		%管存量
Qs = zeros(TimeSteps_Total + 1,1);		%根据时间点处的流量插值计算其它时步点上的流量
Qs(1) = Qs_Opt(1);				%初始点与初始时间点重合
for ii = 1:Secs 					%根据时间点上的值设定整个时间段的流量
	for j = 1:TimeSteps_Per_Sec
		Qs(TimeSteps_Per_Sec*(ii-1)+j+1) = (Qs_Opt(ii+1)-Qs_Opt(ii))*j/TimeSteps_Per_Sec + Qs_Opt(ii);
	end
end
%figure;plot(Qs_Opt);figure;plot(Qs,'r');
Mss = (Den_sta/Area)*Qs;			%起点质量流量密度
Ps = zeros(TimeSteps_Total + 1,1);		%起点压力
Pe = zeros(TimeSteps_Total + 1,1);		%终点压力
Qe = zeros(TimeSteps_Total + 1,1);         	%终点流量
Qe(1) = Ff(1)*Qbasic;
for i = 1:Secs 					%根据时间点上的值设定整个时间段的流量
	for j = 1:TimeSteps_Per_Sec
		Qe(TimeSteps_Per_Sec*(i-1)+j+1) = Qbasic*((Ff(i+1)-Ff(i))*j/TimeSteps_Per_Sec + Ff(i));
	end
end
Mse = Den_sta * Qe/Area;			%终点质量流量密度
%figure;plot(Ff);figure;plot(Qe);

%稳态模拟
tl = Len;               		%管段长度
Ple = 4.5e6;		            	%终点压力
i = SpaceSteps + 1;
Pressure(i) = Ple;		%沿线压力记录
while tl>0              		%稳态模拟
    z = 1 + beta*Ple;   		%压缩因子
    Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(1)^2/C0^2/Din^5;
    Pls = Pls^0.5;
    i = i-1;
    Pressure(i) = Pls;
    tl = tl - dx;
    Ple = Pls;
end
MassFlux = (Den_sta*Qe(1)/Area)*ones(SpaceSteps+1,1);	%构造初始条件-质量流量密度
Pressure = Pressure';		%压力
Storage_Total = 0;		%计算管存量
for ii = 1:SpaceSteps
	Pre_Aver = 2*(Pressure(ii) + Pressure(ii + 1)^2/(Pressure(ii) + Pressure(ii + 1)))/3;
	z = 1 + beta*Pre_Aver;
	Den_Rel = Pre_Aver/z/R/Temp;
	Storage_Sec = Volumn_Sec*Den_Rel;
	Storage_Total = Storage_Total + Storage_Sec;
end
Ps(1) = Pressure(1);
Pe(1) = Pressure(SpaceSteps+1);
Storage(1) = Storage_Total;
%disp(['Initial Total Storage: ' sprintf('%d',Storage_Total)]);
%figure(1);
%plot(Pressure);
%title('Pressure');
%figure(2);
%plot(MassFlux);
%title('MassFlux');

%瞬态模拟
%能耗计算
Total_Com_Consum = 0;
%构造管流控制方程组
creat_transfun_re01(SpaceSteps);
tic;
for i = 2:TimeSteps_Total+1
	tf = @(x)transfun_re01(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Mss(i),Mse(i));	%构造方程
	x0 = zeros(2*SpaceSteps,1);		%准备初值
	x0(1) = Pressure(1);
	for j = 2:SpaceSteps
		x0(2*j-2) = Pressure(j);
		x0(2*j-1) = MassFlux(j);
	end
	x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
	options = optimset('Display','off');
	results = fsolve(tf,x0,options);		%计算
%	results = (1-gama)*x0+gama*results;	%引入阻尼系数，防止计算结果震荡
	Pressure(1) = results(1);		%归档计算结果
	for j = 2: SpaceSteps
		Pressure(j) = results(2*j-2);
		MassFlux(j) = results(2*j-1);
	end
	Pressure(SpaceSteps+1) = results(2*SpaceSteps);
	Ps(i) = results(1);			%管段起点压力
	Pe(i) = results(2*SpaceSteps);		%管段终点压力
	if i < TimeSteps_Total+1
		MassFlux(1) = Mss(i+1);	%引入边界条件
		MassFlux(SpaceSteps+1) = Mse(i+1);
	end
	Storage_Total = 0;			%计算管存量
	for ii = 1:SpaceSteps
		Pre_Aver = 2*(Pressure(ii) + Pressure(ii + 1)^2/(Pressure(ii) + Pressure(ii + 1)))/3;
		z = 1 + beta*Pre_Aver;
		Den_Rel = Pre_Aver/z/R/Temp;
		Storage_Sec = Volumn_Sec*Den_Rel;
		Storage_Total = Storage_Total + Storage_Sec;
	end
	Storage(i) = Storage_Total;
	Quan_Temp = (0.328*Area/Den_sta)*MassFlux(1);	%计算压缩机能耗
	Sec_Com_Consum = dt*Quan_Temp*(2.682*(Pressure(1)/Pin)^0.217 - 2.658);			%该时步压缩机功率
	Total_Com_Consum = Total_Com_Consum + Sec_Com_Consum
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

%存储计算结果
save SimuIm_re03_Results;