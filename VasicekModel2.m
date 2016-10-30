
%%Modelo Vasicek - Riesgo de Mercado UPC %%

function Vasicek=VasicekModel2(B0,B1,T,N,S,move,Var)

% -------------------  Variables ------------------------------------------
%   rt= tasas inicial
%   B0= Beta 0
%   B1= Beta 1
%   N= Cantidad de simulacion (t,T)
%   T= Tiempo de maduración
%   S= Número de simulaciones de tasas
%   Var=Var historico
%   Move= # de intervalos donde el algoritmo toma la tasa real
%   Var_in= Varianza incondicional del modelo Vasicek

% ----------------  Datos a Trabajar  -------------------------------------
% b0	0.000393
% b1	0.991407
% a	0.859%
% b	4.573%
% VAR incondicional	0.17172%
% rt= 5.13%
% VAR	0.00294%
% VasicekModel2(0.000393,0.991407,1,100,100,5,0.0000294)

%--------------------------------------------------------------------------
Data=xlsread('Tes(1-5-10).xls',3);
Contar=length(Data);
T=1;
N=150;
S=100;
move=5;
backtest=Contar-N;
rt=ones(N+1,S);
rt1=ones(N+1,S);


%------- Tasas para cada Move de N resagos---------------------------------

for i=1:move:N+1
rt(i,:)=Data(backtest+i-1,:);
end

for i=1:N+1
rt1(i,:)=Data(backtest+i-1,1);
end
%--------------------------------------------------------------------------

B0=0.000393;
B1=0.991407;
Var=0.00002938;
a=1-B1;
b=B0/(1-B1);
Var_in=(Var)/(1-(B1^2));
dt=T/N;
dW=randn(N+1,S);
unos=ones(N+1,S);
Rate=rt;
RR=rt1;
t= (0:dt:1)';
D= ones(N+1,S);
A= ones(N+1,S);
B= ones(N+1,S);
R= ones(N+1,S);
F= ones(N+1,S);
prom1=ones(N+1,1);


%-------simulación Montecarlo Tasas Cortas -Vasicek Model------------------   

for i=1:N
    for j=1:S
        RR(i+1,j)= RR(i,j)+a*(b-RR(i,j))*dt...
            +sqrt(Var_in)*dW(i,j);
    end
end


for j=1:S
    for z=1:move-1
        for i=z:move:N
   Rate(i+1,j)= Rate(i,j)+a*(b-Rate(i,j))*dt...
       +sqrt(Var_in)*dW(i,j);    
        end    
    end
end
%-- Solucion de la ecuación diferencial con condición final D(t,T)---------

for i=1:N+1
    for j=1:S
        D(i,j)= (1-exp(-a*(T-t(i,1))))/a;
    end
end
%----se eleva al cuadrado cada uno de los elemento de D--------------------

D2=D.^2;

%-- Solucion de la ecuación diferencial con condición de Frontera A(t,T)---

for i=1:N+1
    for j=1:S
        A(i,j) = (1/(a^2))*(D(i,j)-T+t(i,1))*(((a^2)*b)-(0.5*Var_in))...
            -((Var_in*(D2(i,j)))/(4*a));
    end
end

%--------------Precio del Bono---------------------------------------------
for i=1:N+1
    for j=1:S
        B(i,j)=exp(A(i,j)-Rate(i,j)*D(i,j));
    end
end

%-------------Curva Rendimiento Vasicek Model------------------------------

for i=1:N+1
    for j=1:S
        R(i,j)=(-1/(T-t(i,1)))*(log(B(i,j)));
    end
end

%------------- Tasas Forward-----------------------------------------------

for i=1:N+1
    for j=1:S
        F(i,j)= (b-(b-Rate(i,j)))*(exp(-a*(T-t(i,1))))...
            -((Var_in/2)*D2(i,j));
    end
end

%----- Tasas promedio en T-------------------------------------------------

% TasaPromT=mean(Rate(N+1,:));

%----saca el promedio por simulacion---------------------------------------

for i=1:N+1
        prom2(i,1)=mean(Rate(i,:));
end

for i=1:N+1
        prom1(i,1)=mean(RR(i,:));
end



% surf(Rate)
% title('Vasicek Model Short-Term')
% xlabel('Simulations')
% ylabel('Days')
% zlabel('Rate %')

%----Plotters--------------------------------------------------------------
subplot(3,3,2)
plot(0:dt:T,Rate)
title('r(t,T)-Vasicek - Moving')
subplot(3,3,1)
plot(0:dt:T,RR)
title('r(t,T)-Vasicek - Libre')
subplot(3,3,3)
plot(0:dt:T,D(:,1),'b')
title('D(T,T)=0')
subplot(3,3,4)
plot(0:dt:T,A(:,1),'r')
title('A(T,T)=0')
subplot(3,3,5)
plot(0:dt:T,B)
title('B(t,T)-Modelo Vasicek')
% subplot(3,3,6)
% plot(0:dt:T,R)
% title('R(t,T)-Modelo Vasicek')
% subplot(3,3,6)
% plot(0:dt:T,F)
% title('F(t,T)-modelo Vasicek')
subplot(3,3,6)
plot(0:dt:T,rt1(:,1),'b')
title('Tasa Real')
subplot(3,3,7)
plot(0:dt:T,prom1,'r')
title('Simulación 1')
subplot(3,3,8)
plot(0:dt:T,prom2,'black')
title('Simulación 2')
subplot(3,3,9)
plot(0:dt:T,prom1,'r',0:dt:T,rt1(:,1),'b',0:dt:T,prom2,'black')

%----3 graficas superpuestas
% plot(0:dt:T,rt1(:,1),'b',0:dt:T,prom2,'black') %---tasa 1ra y real
% title('Comparación Gráficas')
% legend('Tasa Real','Simulación1','Simulación2')
% subplot(4,3,10)

% Graficos de Box plot 
% boxplot(Rate), hold on
% plot(rt1)
% hold off
% title('Vasicek Model')
% xlabel('Simulations')

% plot(0:dt:T,rt1(:,1),'b',0:dt:T,prom1,'r',0:dt:T,prom2,'black')
% title('Comparación Gráficas - Modelo Vasicek')
% legend('Tasa Real','Simulación 1','Simulación 2')
% xlabel('Days')
% ylabel('Rates')