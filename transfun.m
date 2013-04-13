%状体转移方程

function F = transfun(x,dt,Len,alpha,beta,lamda,Din,Pspre,Pepre,Msspre,Msepre,Psn,Msen)
	F = [(alpha*Psn/(1+beta*Psn)+alpha*x(1)/(1+beta*x(1))-alpha*Pspre/(1+beta*Pspre)-alpha*Pepre/(1+beta*Pepre))/2/dt + (Msen+Msepre-x(2)-Msspre)/2/Len;
	(x(2)+Msen+Msspre+Msepre)/2/dt + ((Msen^2)*(1+beta*x(1))/alpha/x(1)+x(1)+(Msepre^2)*(1+beta*Pepre)/alpha/Pepre+Pepre-(x(2)^2)*(1+beta*Psn)/alpha/Psn-Psn-(Msspre^2)*(1+beta*Pspre)/alpha/Pspre-Pspre)/2/Len + lamda*((x(2)^2)*(1+beta*Psn)/alpha/Psn+(Msen^2)*(1+beta*x(1))/alpha/x(1)+(Msspre^2)*(1+beta*Pspre)/alpha/Pspre+(Msepre^2)*(1+beta*Pepre)/alpha/Pepre)/8/Din];