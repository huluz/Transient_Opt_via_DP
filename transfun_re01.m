%状体转移方程

function F = transfun_re01(x,dt,dx,alpha,beta,lamda,Din,Ppre,Mspre,Mss,Mse)
	F = [(alpha*x(1)/(1+beta*x(1))+alpha*x(2)/(1+beta*x(2))-alpha*Ppre(1)/(1+beta*Ppre(1))-alpha*Ppre(2)/(1+beta*Ppre(2)))/2/dt + (x(3)+Mspre(2)-Mss-Mspre(1))/2/dx;
	(Mss+x(3)-Mspre(1)-Mspre(2))/2/dt + (x(3)^2*(1+beta*x(2))/alpha/x(2)+x(2)+Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2)+Ppre(2)-Mss^2*(1+beta*x(1))/alpha/x(1)-x(1)-Mspre(1)^2*(1+beta*Ppre(1))/alpha/Ppre(1)-Ppre(1))/2/dx + lamda*(Mss^2*(1+beta*x(1))/alpha/x(1)+x(3)^2*(1+beta*x(2))/alpha/x(2)+Mspre(1)^2*(1+beta*Ppre(1))/alpha/Ppre(1)+Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2))/8/Din;
	(alpha*x(2)/(1+beta*x(2))+alpha*x(4)/(1+beta*x(4))-alpha*Ppre(2)/(1+beta*Ppre(2))-alpha*Ppre(3)/(1+beta*Ppre(3)))/2/dt + (x(5)+Mspre(3)-x(3)-Mspre(2))/2/dx;
	(x(3)+x(5)-Mspre(2)-Mspre(3))/2/dt + (x(5)^2*(1+beta*x(4))/alpha/x(4)+x(4)+Mspre(3)^2*(1+beta*Ppre(3))/alpha/Ppre(3)+Ppre(3)-x(3)^2*(1+beta*x(2))/alpha/x(2)-x(2)-Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2)-Ppre(2))/2/dx + lamda*(x(3)^2*(1+beta*x(2))/alpha/x(2)+x(5)^2*(1+beta*x(4))/alpha/x(4)+Mspre(2)^2*(1+beta*Ppre(2))/alpha/Ppre(2)+Mspre(3)^2*(1+beta*Ppre(3))/alpha/Ppre(3))/8/Din;
	(alpha*x(4)/(1+beta*x(4))+alpha*x(6)/(1+beta*x(6))-alpha*Ppre(3)/(1+beta*Ppre(3))-alpha*Ppre(4)/(1+beta*Ppre(4)))/2/dt + (x(7)+Mspre(4)-x(5)-Mspre(3))/2/dx;
	(x(5)+x(7)-Mspre(3)-Mspre(4))/2/dt + (x(7)^2*(1+beta*x(6))/alpha/x(6)+x(6)+Mspre(4)^2*(1+beta*Ppre(4))/alpha/Ppre(4)+Ppre(4)-x(5)^2*(1+beta*x(4))/alpha/x(4)-x(4)-Mspre(3)^2*(1+beta*Ppre(3))/alpha/Ppre(3)-Ppre(3))/2/dx + lamda*(x(5)^2*(1+beta*x(4))/alpha/x(4)+x(7)^2*(1+beta*x(6))/alpha/x(6)+Mspre(3)^2*(1+beta*Ppre(3))/alpha/Ppre(3)+Mspre(4)^2*(1+beta*Ppre(4))/alpha/Ppre(4))/8/Din;
	(alpha*x(6)/(1+beta*x(6))+alpha*x(8)/(1+beta*x(8))-alpha*Ppre(4)/(1+beta*Ppre(4))-alpha*Ppre(5)/(1+beta*Ppre(5)))/2/dt + (x(9)+Mspre(5)-x(7)-Mspre(4))/2/dx;
	(x(7)+x(9)-Mspre(4)-Mspre(5))/2/dt + (x(9)^2*(1+beta*x(8))/alpha/x(8)+x(8)+Mspre(5)^2*(1+beta*Ppre(5))/alpha/Ppre(5)+Ppre(5)-x(7)^2*(1+beta*x(6))/alpha/x(6)-x(6)-Mspre(4)^2*(1+beta*Ppre(4))/alpha/Ppre(4)-Ppre(4))/2/dx + lamda*(x(7)^2*(1+beta*x(6))/alpha/x(6)+x(9)^2*(1+beta*x(8))/alpha/x(8)+Mspre(4)^2*(1+beta*Ppre(4))/alpha/Ppre(4)+Mspre(5)^2*(1+beta*Ppre(5))/alpha/Ppre(5))/8/Din;
	(alpha*x(8)/(1+beta*x(8))+alpha*x(10)/(1+beta*x(10))-alpha*Ppre(5)/(1+beta*Ppre(5))-alpha*Ppre(6)/(1+beta*Ppre(6)))/2/dt + (x(11)+Mspre(6)-x(9)-Mspre(5))/2/dx;
	(x(9)+x(11)-Mspre(5)-Mspre(6))/2/dt + (x(11)^2*(1+beta*x(10))/alpha/x(10)+x(10)+Mspre(6)^2*(1+beta*Ppre(6))/alpha/Ppre(6)+Ppre(6)-x(9)^2*(1+beta*x(8))/alpha/x(8)-x(8)-Mspre(5)^2*(1+beta*Ppre(5))/alpha/Ppre(5)-Ppre(5))/2/dx + lamda*(x(9)^2*(1+beta*x(8))/alpha/x(8)+x(11)^2*(1+beta*x(10))/alpha/x(10)+Mspre(5)^2*(1+beta*Ppre(5))/alpha/Ppre(5)+Mspre(6)^2*(1+beta*Ppre(6))/alpha/Ppre(6))/8/Din;
	(alpha*x(10)/(1+beta*x(10))+alpha*x(12)/(1+beta*x(12))-alpha*Ppre(6)/(1+beta*Ppre(6))-alpha*Ppre(7)/(1+beta*Ppre(7)))/2/dt + (x(13)+Mspre(7)-x(11)-Mspre(6))/2/dx;
	(x(11)+x(13)-Mspre(6)-Mspre(7))/2/dt + (x(13)^2*(1+beta*x(12))/alpha/x(12)+x(12)+Mspre(7)^2*(1+beta*Ppre(7))/alpha/Ppre(7)+Ppre(7)-x(11)^2*(1+beta*x(10))/alpha/x(10)-x(10)-Mspre(6)^2*(1+beta*Ppre(6))/alpha/Ppre(6)-Ppre(6))/2/dx + lamda*(x(11)^2*(1+beta*x(10))/alpha/x(10)+x(13)^2*(1+beta*x(12))/alpha/x(12)+Mspre(6)^2*(1+beta*Ppre(6))/alpha/Ppre(6)+Mspre(7)^2*(1+beta*Ppre(7))/alpha/Ppre(7))/8/Din;
	(alpha*x(12)/(1+beta*x(12))+alpha*x(14)/(1+beta*x(14))-alpha*Ppre(7)/(1+beta*Ppre(7))-alpha*Ppre(8)/(1+beta*Ppre(8)))/2/dt + (x(15)+Mspre(8)-x(13)-Mspre(7))/2/dx;
	(x(13)+x(15)-Mspre(7)-Mspre(8))/2/dt + (x(15)^2*(1+beta*x(14))/alpha/x(14)+x(14)+Mspre(8)^2*(1+beta*Ppre(8))/alpha/Ppre(8)+Ppre(8)-x(13)^2*(1+beta*x(12))/alpha/x(12)-x(12)-Mspre(7)^2*(1+beta*Ppre(7))/alpha/Ppre(7)-Ppre(7))/2/dx + lamda*(x(13)^2*(1+beta*x(12))/alpha/x(12)+x(15)^2*(1+beta*x(14))/alpha/x(14)+Mspre(7)^2*(1+beta*Ppre(7))/alpha/Ppre(7)+Mspre(8)^2*(1+beta*Ppre(8))/alpha/Ppre(8))/8/Din;
	(alpha*x(14)/(1+beta*x(14))+alpha*x(16)/(1+beta*x(16))-alpha*Ppre(8)/(1+beta*Ppre(8))-alpha*Ppre(9)/(1+beta*Ppre(9)))/2/dt + (x(17)+Mspre(9)-x(15)-Mspre(8))/2/dx;
	(x(15)+x(17)-Mspre(8)-Mspre(9))/2/dt + (x(17)^2*(1+beta*x(16))/alpha/x(16)+x(16)+Mspre(9)^2*(1+beta*Ppre(9))/alpha/Ppre(9)+Ppre(9)-x(15)^2*(1+beta*x(14))/alpha/x(14)-x(14)-Mspre(8)^2*(1+beta*Ppre(8))/alpha/Ppre(8)-Ppre(8))/2/dx + lamda*(x(15)^2*(1+beta*x(14))/alpha/x(14)+x(17)^2*(1+beta*x(16))/alpha/x(16)+Mspre(8)^2*(1+beta*Ppre(8))/alpha/Ppre(8)+Mspre(9)^2*(1+beta*Ppre(9))/alpha/Ppre(9))/8/Din;
	(alpha*x(16)/(1+beta*x(16))+alpha*x(18)/(1+beta*x(18))-alpha*Ppre(9)/(1+beta*Ppre(9))-alpha*Ppre(10)/(1+beta*Ppre(10)))/2/dt + (x(19)+Mspre(10)-x(17)-Mspre(9))/2/dx;
	(x(17)+x(19)-Mspre(9)-Mspre(10))/2/dt + (x(19)^2*(1+beta*x(18))/alpha/x(18)+x(18)+Mspre(10)^2*(1+beta*Ppre(10))/alpha/Ppre(10)+Ppre(10)-x(17)^2*(1+beta*x(16))/alpha/x(16)-x(16)-Mspre(9)^2*(1+beta*Ppre(9))/alpha/Ppre(9)-Ppre(9))/2/dx + lamda*(x(17)^2*(1+beta*x(16))/alpha/x(16)+x(19)^2*(1+beta*x(18))/alpha/x(18)+Mspre(9)^2*(1+beta*Ppre(9))/alpha/Ppre(9)+Mspre(10)^2*(1+beta*Ppre(10))/alpha/Ppre(10))/8/Din;
	(alpha*x(18)/(1+beta*x(18))+alpha*x(20)/(1+beta*x(20))-alpha*Ppre(10)/(1+beta*Ppre(10))-alpha*Ppre(11)/(1+beta*Ppre(11)))/2/dt + (Mse+Mspre(11)-x(19)-Mspre(10))/2/dx;
	(x(19)+Mse-Mspre(10)-Mspre(11))/2/dt + (Mse^2*(1+beta*x(20))/alpha/x(20)+x(20)+Mspre(11)^2*(1+beta*Ppre(11))/alpha/Ppre(11)+Ppre(11)-x(19)^2*(1+beta*x(18))/alpha/x(18)-x(18)-Mspre(10)^2*(1+beta*Ppre(10))/alpha/Ppre(10)-Ppre(10))/2/dx + lamda*(x(19)^2*(1+beta*x(18))/alpha/x(18)+Mse^2*(1+beta*x(20))/alpha/x(20)+Mspre(10)^2*(1+beta*Ppre(10))/alpha/Ppre(10)+Mspre(11)^2*(1+beta*Ppre(11))/alpha/Ppre(11))/8/Din;
	];