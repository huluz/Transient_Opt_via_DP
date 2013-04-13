%采用隐式算法进行非稳态模拟

clc;clear;
mol = 17.44;            %气体分子质量
R = 8314;               %通用气体常数
PreC = 46.59E5;         %临界压力
TempC = 198.19;         %临界温度
Temp = 15+273.15;       %输气温度
Temp_sta = 20+273.15;   %标况温度
Pre_sta = 101325;       %标况压力
Den_air = 1.209;        %空气密度
Rel_Den = 0.6;          %天然气相对密度
Den_sta = Rel_Den*Den_air;  %天然气标况密度
alpha = Den_sta*Temp_sta/Pre_sta/Temp;  %通用气体方程参数
beta = 0.257/PreC - 0.533*TempC/PreC/Temp;  %AGA气体状态方程参数

Len = 100E3;            %管段长度
Time = 24*3600;			%模拟时长
dt = 3600;              %时步
dx = 2E3;               %稳态模拟空间步长
lamda = 0.008;          %摩阻系数
Din = 0.6096;           %管段内径

TimeSteps = Time/dt; 	%时步数
Ps = 7.3E6*ones(24,1);	%管段起点压力
Qs = zeros(24,1);		%起点压力
Mss = zeros(24,1);		%起点流量
Pe = zeros(24,1);		%终点压力
Qbasic = 65;            %流量基数
Ff = [0.2;0.15;0.1;0.25;0.35;0.58;1.2;1.3;1;0.97;0.85;1.65;2;1;0.8;0.65;...
    1.15;1.9;2.8;2.2;1.2;0.85;0.5;0.35];    %小时流量不均匀系数
Qe = Qbasic*Ff;         %终点流量
Mse = Den_sta*Qe;		%终点质量流量密度

C0 = 0.03848;           %稳态模拟公式系数

%稳态模拟
tl = Len;               %管段长度
Pls = Ps(1);            %起点压力
while tl>0              %稳态模拟
    z = 1 + beta*Pls;   %压缩因子
    Ple = Pls^2 - lamda*z*Rel_Den*Temp*dx*Qe(24)^2/C0^2/Din^5;
    Ple = Ple^0.5;
    tl = tl - dx;
    Pls = Ple;
end
Ps_ini = Ps(1);            	%起点初始压力
Mss_ini = Den_sta*Qe(24);	%初始质量流量密度
Pe_ini = Ple;       		%终点质量流量密度
Mse_ini = Den_sta*Qe(24);  	%初始质量流量密度

%瞬态模拟
%第一时步
tf = @(x)transfun(x,dt,Len,alpha,beta,lamda,Din,Ps_ini,Pe_ini,Mss_ini,Mse_ini,Ps(1),Mse(1));
x0 = [Pe_ini,Mss_ini];
%options = optimset('MaxFunEvals',1E4,'MaxIter',1E4);
results = fsolve(tf,x0);
x=1;

for i = 2:TimeSteps
	tf = @(Mssn,Pen) tranfun(dt,Len,alpha,beta,lamda,Din,Pspre,Pepre,Msspre,Msepre,Psn,Pen,Mssn,Msen);
end