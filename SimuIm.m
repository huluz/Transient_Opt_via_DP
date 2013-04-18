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
Area = 0.25*pi*Din^2;		%管段横截面积
C0 = 0.03848;           		%稳态模拟公式系数
Time = 24*3600;		%模拟时长
Time_Sec = 3600;		%单个时间段时长
dt = 60;              		%时步
TimeSteps_Total = Time / dt; 	%总时步数
TimeSteps_Per_Sec = Time_Sec / dt;		%每个时间段的步长数
dx = 10E3;               		%空间步长
SpaceSteps = Len / dx;		%空间分段数
gama = 0.95;			%阻尼系数
tf_str = '';			%状态转移方程字符串

%边界条件
Ps = 7.3E6*ones(TimeSteps_Total,1);		%管段起点压力
Qs = zeros(TimeSteps_Total,1);		%起点流量
Mss = zeros(TimeSteps_Total,1);		%起点质量流量密度
Pe = zeros(TimeSteps_Total,1);		%终点压力
Qbasic = 65;           		 		%流量基数
Ff = [0.2;0.15;0.1;0.25;0.35;0.58;1.2;1.3;1;0.97;0.85;1.65;2;1;0.8;0.65;...
    1.15;1.9;2.8;2.2;1.2;0.85;0.5;0.35];    		%小时流量不均匀系数
Qe = zeros(TimeSteps_Total,1);         		%终点流量
for i = 1:Time/Time_Sec 			%根据时间点上的值设定整个时间段的流量
	Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
end
Mse = Den_sta * Qe/Area;			%终点质量流量密度

%稳态模拟
tl = Len;               		%管段长度
Pls = Ps(1);            		%起点压力
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
%figure(1);
%plot(Pressure);
%title('Pressure');
%figure(2);
%plot(MassFlux);
%title('MassFlux');

%瞬态模拟
%构造管流控制方程组
fid = fopen('transfun.m');	%打开状态转移方程文件
fprintf(fid,'%s','function F = transfun(x,dt,dx,alpha,beta,lamda,Din,Ppre,Mspre,Ps,Mse)\n');
fprintf(fid,'%s','F = [(alpha*Ps/(1+beta*Ps)+alpha*x(2)/(1+beta*x(2))-alpha*Ppre(1)/(1+beta*Ppre(1))-alpha*Ppre(2)/(1+beta*Ppre(2)))/2/dt + (x(3)+Mspre(2)-x(1)-Mspre(1))/2/dx;\n');
fprintf(fid,'%s','(x(1)+x(3)-Mspre(1)-Mspre(2))/2/dt + (x(3)^2*(1+beta*x(2))/alpha/x(2)+x(2)+Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2)+Ppre(2)-x(1)^2*(1+beta*Ps)/alpha/Ps-Ps-Mspre(1)^2*(1+beta*Ppre(1))/alpha/Ppre(1)-Ppre(1))/2/dx + lamda*(x(1)^2*(1+beta*Ps)/alpha/Ps+x(3)^2*(1+beta*x(2))/alpha/x(2)+Mspre(1)^2*(1+beta*Ppre(1))/alpha/Ppre(1)+Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2))/8/Din;\n');
for i = 2:SpaceSteps-1
	fprintf(fid,'%s','(alpha*x(2*i)/(1+beta*x(2*i-2))+alpha*x(2*i)/(1+beta*x(2*i))-alpha*Pressure(i)/(1+beta*Pressure(i))-alpha*Pressure(i+1)/(1+beta*Pressure(i+1)))/2/dt + (x(2*i+1)+MassFlux(i+1)-x(2*i-1)-MassFlux(i))/2/dx;');
	fprintf(fid,'%s','(x(2*i-1)+x(2*i+1)-MassFlux(i)-MassFlux(i+1))/2/dt + ((x(2*i+1)^2)*(1+beta*x(2*i-2))/alpha/x(2*i-2)+x(2*i)+(MassFlux(i+1)^2)*(1+beta*Pressure(i+1))/alpha/Pressure(i+1)+Pressure(i+1)-(x(2*i-1)^2)*(1+beta*x(2*i-2))/alpha/x(2*i-2)-x(2*i-2)-(MassFlux(i)^2)*(1+beta*Pressure(i))/alpha/Pressure(i)-Pressure(i))/2/dx + lamda*((x(2*i+1)^2)*(1+beta*x(2*i-2))/alpha/x(2*i-2)+(MassFlux(i+1)^2)*(1+beta*Pressure(i+1))/alpha/Pressure(i+1)+(x(2*i-1)^2)*(1+beta*x(2*i-2))/alpha/x(2*i-2)+(MassFlux(i)^2)*(1+beta*Pressure(i))/alpha/Pressure(i))/8/Din')
end
fprintf(fid,'%s','(alpha*x(2*SpaceSteps-2)/(1+beta*x(2*SpaceSteps-2))+alpha*x(2*SpaceSteps)/(1+beta*x(2*SpaceSteps))-alpha*Pressure(SpaceSteps)/(1+beta*Pressure(SpaceSteps))-alpha*Pressure(SpaceSteps+1)/(1+Pressure(SpaceSteps+1)))/2/dt + (Mse(1)+MassFlux(SpaceSteps+1)-x(2*SpaceSteps-1)-MassFlux(SpaceSteps))/2/dx;');
fprintf(fid,'%s','x(2*SpaceSteps-1)+Mse(1)-MassFlux(SpaceSteps)-MassFlux(SpaceSteps+1))/2/dt + ((Mse(1)^2)*(1+beta*x(2*SpaceSteps))/alpha/x(2*SpaceSteps)+x(2*SpaceSteps)+(MassFlux(SpaceSteps+1)^2)*(1+beta*Pressure(SpaceSteps+1))/alpha/Pressure(SpaceSteps+1)+Pressure(SpaceSteps+1)-(x(2*i-1)^2)*(1+beta*x(2*SpaceSteps-2))/alpha/x(2*SpaceSteps-2)-x(2*SpaceSteps-2)-(MassFlux(SpaceSteps)^2)*(1+beta*Pressure(SpaceSteps))/alpha/Pressure(SpaceSteps)-Pressure(SpaceSteps-1))/2/dx + lamda*((Mse(1)^2)*(1+beta*x(2*SpaceSteps))/alpha/x(2*SpaceSteps)+(MassFlux(SpaceSteps+1)^2)*(1+beta*Pressure(SpaceSteps+1))/alpha/Pressure(SpaceSteps+1)+(x(2*i-1)^2)*(1+beta*x(2*SpaceSteps-2))/alpha/x(2*SpaceSteps-2)+(MassFlux(SpaceSteps)^2)*(1+beta*Pressure(SpaceSteps))/alpha/Pressure(SpaceSteps))/8/Din];');
tic;
for i = 1:TimeSteps_Total
	tf = @(x)transfun(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Ps(i),Mse(i));	%构造方程
	x0 = zeros(2*SpaceSteps,1);	%准备初值
	x0(1) = MassFlux(1);
	for j = 2:SpaceSteps
		x0(2*j-2) = Pressure(j);
		x0(2*j-1) = MassFlux(j);
	end
	x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
	%options = optimset('MaxFunEvals',1E4,'MaxIter',1E4,'Display','iter-detailed');
	options = optimset('Display','off');
	results = fsolve(tf,x0,options);		%计算
	results = (1-gama)*x0+gama*results;	%引入阻尼系数，防止计算结果震荡
	MassFlux(1) = results(1);	%归档计算结果
	for j = 2: SpaceSteps
		Pressure(j) = results(2*j-2);
		MassFlux(j) = results(2*j-1);
	end

	Pressure(SpaceSteps+1) = results(2*SpaceSteps);
	Mss(i) = results(1);		%管段起点质量流量
	Pe(i) = results(2*SpaceSteps);	%管段终点压力
	if i < TimeSteps_Total
		Pressure(1) = Ps(i+1);		%引入边界条件
		MassFlux(SpaceSteps+1) = Mse(i+1);
	end
	if mod(i,10) == 0
		%计算过程图形化
		figure(3);			%终点压力变化
		plot(Pe(1:i));
		title('Pressure at the Outlet');
		Qs(1:i) = (Area/Den_sta)*Mss(1:i);		%起点流量变化
		figure(4);
		xaxes = 1:i;
		plot(xaxes,Qs(1:i),xaxes,Qe(1:i));
		title('Quantity');
		legend('Quantity at the Inlet','Quantity at the Outlet');
	end
end
toc;

%计算结果可视化
figure(3);			%终点压力变化
plot(Pe);
title('Pressure at the Outlet');
Qs = (Area/Den_sta)*Mss;	%起点流量变化
figure(4);
xaxes = 1:TimeSteps_Total;
plot(xaxes,Qs,xaxes,Qe);
title('Quantity at the Inlet');
legend('Quantity at the Inlet','Quantity at the Outlet');