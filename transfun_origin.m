%状体转移方程

function F = transfun_origin(x,dt,Len,alpha,beta,lamda,Din,Pspre,Pepre,Msspre,Msepre,Psn,Msen)
	F = [(alpha*Psn/(1+beta*Psn)+alpha*Pen/(1+beta*Pen)-alpha*Pspre/(1+beta*Pspre)-alpha*Pepre/(1+beta*Pepre))/2/dt + (Msen+Msepre-Mssn-Msspre)/2/Len;
	(Mssn+Msen-Msspre-Msepre)/2/dt + ((Msen^2)*(1+beta*Pen)/alpha/Pen+Pen+(Msepre^2)*(1+beta*Pepre)/alpha/Pepre+Pepre-(Mssn^2)*(1+beta*Psn)/alpha/Psn-Psn-(Msspre^2)*(1+beta*Pspre)/alpha/Pspre-Pspre)/2/Len + lamda*((Mssn^2)*(1+beta*Psn)/alpha/Psn+(Msen^2)*(1+beta*Pen)/alpha/Pen+(Msspre^2)*(1+beta*Pspre)/alpha/Pspre+(Msepre^2)*(1+beta*Pepre)/alpha/Pepre)/8/Din];