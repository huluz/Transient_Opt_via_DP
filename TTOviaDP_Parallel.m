% ����DP�㷨��������ܵ�ĩ�η���̬�Ż�����

% �����ʼ��
clc;clear;

%������ز���
mol = 17.44;			%��Է�������
R = 8314;			%ͨ�����峣��
PreC = 46.59E5;		%�������ٽ�ѹ��
TempC = 198.19;		%�������ٽ��¶�
Temp = 15 + 273.15;		%�����¶�
Den_air = 1.209;		%����¿����ܶȣ�kg/m^3
Temp_sta = 20 + 273.15;	%����¶�
Pre_sta = 101325;		%���ѹ��
Rel_Den = 0.6;			%����ܶ�
Den_sta = Rel_Den * Den_air;	%��������ܶ�
alpha = Den_sta*Temp_sta/Pre_sta/Temp;	%ͨ������״̬���̲���
beta = 0.257/PreC - 0.533*TempC/PreC/Temp;	%AGA����״̬���̲���

%�ܶβ���
Len = 100E3;			%�ܵ����ȣ�m
Pin = 3.5E6;			%ѹ�������ѹ��
lamda = 0.008;			%Ħ��ϵ��
Din = 0.6096;			%���ھ���m
Area = 0.25*pi*Din^2;		%�ܶκ�����
Pe_min = 4.5E6;		%�ܶγ����������ѹ��
Ps_max = 7.3E6;		%�ܶν����������ѹ��
Qbasic = 65;			%����������m^3/s
Ff = [0.2; 0.15; 0.1; 0.25; 0.35; 0.58; 1.2; 1.3; 1; 0.97; 0.85; 1.65; 2; 1; 0.8; 0.65; ...
	1.15; 1.9; 2.8; 2.2; 1.2; 0.85; 0.5; 0.35];	%Сʱ����������ϵ��

%�������
Time = 24*3600;		%�Ż�ʱ��Σ�s
Time_Per_Sec = 3600;		%�����߽������趨ʱ����s
dt = 1800;			%ʱ����s
Time_Secs = Time / Time_Per_Sec;	%ʱ�����
TimeSteps_Total = Time / dt;	%��ʱ����
TimeSteps_Per_Sec = Time_Per_Sec / dt;	%ÿ��ʱ��ΰ�����ʱ����
dp = 0.2E6;			%״̬�����ռ���ɢ����
dx = 10E3;			%�ռ䲽����m
SpaceSteps = Len / dx;		%�ռ�ֶ���
creat_transfun(SpaceSteps);	%����״̬ת�Ʒ���

%ģ�����
C0 = 0.03848;			%��̬ģ�⹫ʽ����
a = [4.4889;0.071;-0.009];	%ѹ����ѹͷ-��������ϵ��
b = [9.0719;2.2055;-0.0161];	%ѹ����Ч��-��������ϵ��
Pressure_ini = zeros(SpaceSteps+1,1);		%����ѹ���ֲ�
MassFlux_ini = zeros(SpaceSteps+1,1);	%�������������ֲ�

%�Ż��������

%�߽�����
Qe = zeros(TimeSteps_Total,1);		%�յ�����
for i = 1:Time_Secs			%����ʱ����ϵ�ֵ�趨����ʱ��ε�����
	%Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
	for j = 1:TimeSteps_Per_Sec
		if i == 1
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(Time_Secs))*j/TimeSteps_Per_Sec + Ff(Time_Secs));
		else
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(i-1))*j/TimeSteps_Per_Sec + Ff(i-1));
		end
	end
end
Mse = Den_sta * Qe/Area;	%�յ����������ܶ�

%DP�㷨��ʼ��-���ó�ʼ״̬
%��̬ģ��
tl = Len;			%�ܶγ���
Ple = Pe_min;			%���ѹ��
Pressure_ini(SpaceSteps+1) = Ple;%����ѹ����¼
i = SpaceSteps;
while tl>0 			%��̬ģ��
	z = 1 + beta*Ple;	%ѹ������
	Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(TimeSteps_Total)^2/C0^2/Din^5;
	Pls = Pls^0.5;
	Pressure_ini(i) = Pls;
	i = i - 1;
	tl = tl - dx;
	Ple = Pls;
end
MassFlux_ini = (Den_sta*Qe(TimeSteps_Total)/Area)*ones(SpaceSteps+1,1);	%�����ʼ����-���������ܶ�
Pressure_ini = Pressure_ini';		%ѹ��
%figure(1);
%plot(Pressure_ini);
%title('Pressure_ini');
%figure(2);
%plot(MassFlux_ini);
%title('MassFlux_ini');

%DP�㷨����Ż�����

%����DP�㷨��ʼ״̬
Results_Now = [Pressure_ini,MassFlux_ini',zeros(1,Time_Secs),0];	%�趨��ʼ״̬
%Storage = 0;								%�ܴ���
Volumn_Sec = Area*dx;
%for ii = 1:SpaceSteps
%	Pre_Aver = (Pressure_ini(ii) + Pressure_ini(ii+1))/2;
%	z = 1 + beta*Pre_Aver;
%	Den_Rel = Pre_Aver/z/R/Temp;
%	Storage_Sec = Volumn_Sec*Den_Rel;
%	Storage = Storage + Storage_Sec;
%end
Col_Num_Per_Rec = 2*(SpaceSteps + 1) + Time_Secs + 1;		%ÿ����¼���������
Ms_Rec_Start_Num = SpaceSteps + 2;					%����������¼��ʼλ��
Desicision_Rec_Start_Num = 2*(SpaceSteps + 1) + 1;			%�������м�¼��ʼ����

%˳����ƹ���
for i = 1:5
	startmatlabpool;
	disp('=========================================');
	disp(['Time: ' sprintf('%d',i)]);
	tic;				%��ʱ
	%ȷ�����߱�����Χ
	tl = Len;			%�ܶγ���
	Ple = Pe_min;			%���ѹ��
	while tl>0 			%��̬ģ��
		z = 1 + beta*Ple;	%ѹ������
		Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(i*TimeSteps_Per_Sec)^2/C0^2/Din^5;
		Pls = Pls^0.5;
		tl = tl - dx;
		Ple = Pls;
	end
	Ps_min = Pls;

	%��ɢ���߱���������
	Ps_avai_Num = 0; 		%������ɢ��õ��ľ��߱�����ȡֵ��Ŀ
	while Ps_max - Ps_avai_Num*dp > Ps_min
		Ps_avai_Num = Ps_avai_Num + 1;
	end
	if Ps_avai_Num <= 0
		sprintf('%s\n','Desicision is Empty!');
		return;
	end
	%Ps_avai_Num
	Ps_avai = zeros(Ps_avai_Num,1);	%������߱�����ȡֵ
	for j = 1:Ps_avai_Num
		Ps_avai(j) = Ps_max - (j - 1)*dp;
	end

	%ת����һ�׶εļ���������׼����������������ڴ�ռ�
	Results_Pre = Results_Now;			%ת����һ�׶εļ�����
%	Results_Pre_Pressure = Results_Pre(:, 1:Ms_Rec_Start_Num - 1);				%��ּ�¼��ѹ����¼
%	Results_Pre_MassFlux = Results_Pre(:, Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1);	%����������¼
%	Results_Pre_Decision = Results_Pre(:, Desicision_Rec_Start_Num:Col_Num_Per_Rec - 1);	%���߱�����¼
%	Results_Pre_Com_Consum = Results_Pre(:, Col_Num_Per_Rec);				%ѹ�����ܺļ�¼
	States_Pre_Num = size(Results_Pre,1);	%���׶�/��һʱ��״̬��
	States_Now_Num_Temp = States_Pre_Num*Ps_avai_Num;			%��һ�׶�/��ʱ��״̬����������
	Results_Now_Temp = zeros(States_Now_Num_Temp,Col_Num_Per_Rec);	%׼�����������ڴ�ռ�
	Results_Now_Temp_Pressure = zeros(States_Now_Num_Temp, Ms_Rec_Start_Num - 1);				%ѹ����¼
	Results_Now_Temp_MassFlux = zeros(States_Now_Num_Temp, Desicision_Rec_Start_Num - Ms_Rec_Start_Num);	%����������¼
	Results_Now_Temp_Decision = zeros(States_Now_Num_Temp, Col_Num_Per_Rec - Desicision_Rec_Start_Num);	%���߹��̼�¼
	Results_Now_Temp_Decision_TT = zeros(States_Now_Num_Temp,1);
	Results_Now_Temp_Com_Consum = zeros(States_Now_Num_Temp, 1);						%ѹ�����ܺļ�¼
	Storage = zeros(States_Now_Num_Temp,1);	%�ܴ���

	%����״̬����߱�����ϣ�������ģ��
	parfor Rec_Num = 1:States_Now_Num_Temp 		%���չ���״̬������ѭ��
		%����ѭ���������õ�״̬����߱�������
		n = mod(Rec_Num, Ps_avai_Num) ;	%���߱�������
		if n == 0 
			n = Ps_avai_Num;
			m = Rec_Num/Ps_avai_Num;
		else
			m = (Rec_Num - n)/Ps_avai_Num + 1;	%״̬����
		end

		%���ƾ��߹���
		Results_Now_Temp_Decision(Rec_Num, :) = Results_Pre(m, Desicision_Rec_Start_Num:Col_Num_Per_Rec-1);
		Results_Now_Temp_Decision_TT(Rec_Num) = Ps_avai(n);			%�趨���߱���ֵ
		Pressure = Results_Pre(m,1:Ms_Rec_Start_Num - 1);					%�趨ѹ����ʼֵ��׼��ģ��
		MassFlux = Results_Pre(m,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1);	%�趨����������ʼֵ

		%ģ�����
		Ps_Pre = Results_Pre(m,1);	%ȡ����һʱ���ܶ���ڴ�ѹ��
		Com_Consum = 0;		%ѹ�����ܺ�
		%����ʱ����ڿ��ܴ��ڶ��ʱ������Ҫ���ж��ģ��
		for l = 1:TimeSteps_Per_Sec
			Ps_Simu = (Ps_avai(n) - Ps_Pre)*l/TimeSteps_Per_Sec + Ps_Pre;	%ͨ����ֵ���������߽�����
			tf = @(x)transfun(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Ps_Simu,Mse((i-1)*TimeSteps_Per_Sec+l));	%���췽��
			x0 = zeros(2*SpaceSteps,1);	%׼����ֵ
			x0(1) = MassFlux(1);
			for j = 2:SpaceSteps
				x0(2*j-2) = Pressure(j);
				x0(2*j-1) = MassFlux(j);
			end
			x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
			options = optimset('Display','off');
			results = fsolve(tf,x0,options);	%����

			MassFlux(1) = results(1);	%�鵵������
			for j = 2:SpaceSteps
				Pressure(j) = results(2*j-2);
				MassFlux(j) = results(2*j-1);
			end
			Pressure(SpaceSteps+1) = results(2*SpaceSteps);
			Pressure(1) = Ps_Simu;
			MassFlux(SpaceSteps+1) = Mse((i-1)*TimeSteps_Per_Sec+l);

			%����ѹ�����ܺ�
			if Pressure(SpaceSteps+1) < Pe_min || min(MassFlux) < 0
				Results_Now_Temp_Com_Consum(Rec_Num) = -1;
				Com_Consum = 0;
			else
				Com_Consum = Com_Consum + dt*MassFlux(1)*Area*(Pressure(1)/Pin)^0.8;
			end
		end
		%�鵵������
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

	%���������
	for Rec_Num = 1:States_Now_Num_Temp
		Results_Now_Temp(Rec_Num,1:Ms_Rec_Start_Num - 1) = Results_Now_Temp_Pressure(Rec_Num,:);				%ѹ����¼
		Results_Now_Temp(Rec_Num,Ms_Rec_Start_Num:Desicision_Rec_Start_Num - 1) = Results_Now_Temp_MassFlux(Rec_Num,:);	%����������¼
		Results_Now_Temp_Decision(Rec_Num,i) = Results_Now_Temp_Decision_TT(Rec_Num);						%���߹��̼�¼
		Results_Now_Temp(Rec_Num,Desicision_Rec_Start_Num:Col_Num_Per_Rec-1) = Results_Now_Temp_Decision(Rec_Num,:);
		Results_Now_Temp(Rec_Num,Col_Num_Per_Rec) = Results_Now_Temp_Com_Consum(Rec_Num);					%�ܺļ�¼
	end

	%���ݹܴ�����״̬�ռ�
	Red_States = 0;	%����״̬��
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

	%�޳������з���
	Bad_Recs = 0;		%ͳ�ƿ��з�����Ŀ
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
	l = 1;m=1;		%�鵵���з���
	Results_Now = zeros(Good_Recs,Col_Num_Per_Rec);
	while l <= Good_Recs
		if Results_Now_Temp(m,Col_Num_Per_Rec) ~= -1
			Results_Now(l,:) = Results_Now_Temp(m,:);
			l = l + 1;
		end
		m = m + 1;
	end

	%������̿��ӻ�
	disp(['Total Combinations: ' sprintf('%d',States_Now_Num_Temp)]);
	disp(['Good Recodes: ' sprintf('%d',Good_Recs)]);
	disp(['Bad Recodes; ' sprintf('%d',Bad_Recs)]);
	disp(['Redundancy Recodes: ' sprintf('%d', Red_States)]);
	toc;			%����ʱ��μ���ʱ��
	disp('=========================================');
	disp('');
end

%ɸѡ���ŷ���
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

%ͼ�λ�������
%figure(3);
%plot(Opt_Des);
%title('Optimum Control Strategy');
%legend('Volumn quantity at the inlet');