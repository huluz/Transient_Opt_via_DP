% 验证算法�?��用的状�?转移方程

%初始�?
clc;clear;
mol = 17.44;            %相对分子质量
R = 8314;               %通用气体常数
PreC = 46.59E5;         %气体视临界压�?
TempC = 198.19;         %气体视临界温�?
Temp = 15+273.15;       %输气温度
Temp_sta = 20+273.15;   %标况温度
Pre_sta = 101325;       %标况压力
Den_air = 1.209;        %空气标况下的密度，kg/m^3
Rel_Den = 0.6;          %相对密度
Den_sta = Rel_Den*Den_air;  %标况气体密度
alpha = Den_sta*Temp_sta/Pre_sta/Temp;  %通用气体状�?方程参数
beta = 0.257/PreC - 0.533*TempC/PreC/Temp;  %AGA气体状�?方程参数

Ps = 7.3E6*ones(24,1);	%管段进口压力
Qs = zeros(24,1);		%管段起点流量-待求
Mss = zeros(24,1);		%管段起点质量流量-待求
Pe = zeros(24,1);		%管段终点压力-待求
Qbasic = 65;            %基础流量，m^3/s
Ff = [0.2;0.15;0.1;0.25;0.35;0.58;1.2;1.3;1;0.97;0.85;1.65;2;1;0.8;0.65;...
    1.15;1.9;2.8;2.2;1.2;0.85;0.5;0.35];    %小时不均�?��气系�?
Qe = Qbasic*Ff;         %管段末端流量�?��
Mse = Den_sta*Qe;		%管段起点质量流量密度

Len = 100E3;            %管道长度，m
dt = 3600;              %时间步长，s
dx = 2E3;               %
dp = 0.1E6;             %状�?变量空间离散步长
lamda = 0.008;          %摩阻系数
Din = 0.6096;           %管内径，m
a = [4.4889;0.071;-0.009];  %压缩机压�?流量曲线系数
b = [9.0719;2.2055;-0.0161];%压缩机效�?流量曲线系数

C0 = 0.03848;           %稳�?模拟公式参数

%初始状�?
tl = Len;               %稳�?模拟
Pls = Ps(1);            %稳�?模拟划分管段末端压力
while tl>0              %求解压气站出口压�?
    z = 1 + beta*Pls;   %管段终点压缩因子
    Ple = Pls^2 - lamda*z*Rel_Den*Temp*dx*Qe(24)^2/C0^2/Din^5;
    Ple = Ple^0.5;
    tl = tl - dx;
    Pls = Ple;
end
Ps_ini = Ps(1);            	%管段起点压力
Mss_ini = Den_sta*Qe(24);	%质量流量密度
Pe_ini = Ple;       		%管段终点压力
Mse_ini = Den_sta*Qe(24);  	%质量流量密度

%�?��瞬�?模拟
Cu = 2*dt*((Mss_ini^2)*(1+beta*Ps_ini)/alpha/Ps_ini+Ps_ini-(Mse_ini^2)*(1+beta*Pe_ini)/alpha/Pe_ini-Pe_ini)/Len;
Du = dt*(lamda*(Mss_ini^2)*(1+beta*Ps_ini)/alpha/Ps_ini+lamda*(Mse_ini^2)*(1+beta*Pe_ini)/alpha/Pe_ini)/2/Din;
Mss(1) = Mss_ini + Mse_ini - Mse(1) + Cu - Du;
Au = alpha*Ps_ini/(1+beta*Ps_ini) + alpha*Pe_ini/(1+beta*Pe_ini) - alpha*Ps(1)/(1+beta*Ps(1)) + 2*dt*(Mss_ini - Mse_ini)/Len;
Pe(1) = Au/(alpha - Au*beta);

for i = 2:24			
	Cu = 2*dt*((Mss(i-1)^2)*(1+beta*Ps(i-1))/alpha/Ps(i-1) + Ps(i-1) - (Mse(i-1)^2)*(1+beta*Pe(i-1))/alpha/Pe(i-1) - Pe(i-1));
	Du = dt*(lamda*(Mss(i-1)^2)*(1+beta*Ps(i-1))/alpha/Ps(i-1) + lamda*(Mse(i-1)^2)*(1+beta*Pe(i-1))/alpha/Pe(i-1));
	Mss(i) = Mss(i-1) + Mse(i-1) - Mse(i) + Cu - Du;
	Au = alpha*Ps(i-1)/(1+beta*Ps(i-1)) + alpha*Pe(i-1)/(1+beta*Pe(i-1)) - alpha*Ps(i)/(1+beta*Ps(i)) + 2*dt*(Mss(i-1)-Mse(i-1))/Len;
	Pe(i) = Au/(alpha-Au*beta);
end

%作图
figure(1);
plot(Pe);
figure(2);
plot(Mss/Den_sta);