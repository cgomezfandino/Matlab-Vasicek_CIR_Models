
%%Modelo CIR con el MGM - Riesgo de Mercado UPC %%

function CIR=CIRModel2(B0,B1,T,N,S,move,Var)

% -------------------  Variables ------------------------------------------
%   rt= tasas inicial
%   B0= Beta 0
%   B1= Beta 1
%   N= Cantidad de simulacion (t,T)
%   T= Tiempo de maduración
%   S= Número de simulaciones de tasas
%   Var=Var historico
%   Move= # de intervalos donde el algoritmo toma la tasa real

% ----------------  Datos a Trabajar  -------------------------------------
% b0	0.000393
% b1	0.991407
% a	0.859%
% b	4.573%
% VAR incondicional	0.17172%
% rt= 5.13%
% VAR	0.00294%
% CIRModel2(0.000393,0.991407,1,100,100,5,0.0000294)

%--------------------------------------------------------------------------
Data=xlsread('Tes(1-5-10).xls',3);
Contar=length(Data);
% T=1;
% N=150;
% S=100;
% move=5;
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

% B0=0.000393;
% B1=0.991407;
% Var=0.00002938;
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
gamma=0.5*sqrt((a^2)+(2*Var));


%------------------simulaciones--------------------------------------------


%-------simulación Montecarlo Tasas Cortas -CIR Model----------------------   

for i=1:N
    for j=1:S
        RR(i+1,j)= RR(i,j)+a*(b-RR(i,j))*dt...
            +sqrt(Var)*sqrt(RR(i,j))*dW(i,j);
    end
end


for j=1:S
    for z=1:move-1
        for i=z:move:N
   Rate(i+1,j)= Rate(i,j)+a*(b-Rate(i,j))*dt...
       +sqrt(Var)*sqrt(Rate(i,j))*dW(i,j);    
        end    
    end
end


%------- senh--------------------------------------------------------------
for i=1:N+1
    for j=1:S
        senh(i,j)=(exp(gamma*(T-t(i,1)))-...
            exp(-gamma*(T-t(i,1))))/2;
    end
end

%------cosh----------------------------------------------------------------
for i=1:N+1
    for j=1:S
        cosh(i,j)=(exp(gamma*(T-t(i,1)))+...
            exp(-gamma*(T-t(i,1))))/2;
    end
end


%-- Solucion de la ecuación diferencial con condición final D(t,T)---------

for i=1:N+1
    for j=1:S
       D(i,j)=senh(i,j)/((gamma*cosh(i,j))+(0.5*a*senh(i,j)));
    end
end

%------ Otra forma de hallar D(t,T)----------------------------------------
% for i=1:N+1
%     for j=1:S
%         D1(i,j)=(2*(exp(sqrt((a^2)+(2*Var_in))*(T-t(i,1)))-1))...
%             /((a+(sqrt(a^2+2*Var_in)))*(exp(sqrt((a^2)+(2*Var_in))...
%             *(T-t(i,1)))-1)+(2*sqrt((a^2)+(2*Var_in))));
%     end
% end

%-- Solucion de la ecuación diferencial con condición de Frontera A(t,T)---

for i=1:N+1
    for j=1:S
        A(i,j)=(2*a*b/Var)*log((gamma*exp(0.5*a*(T-t(i,1))))...
            /((gamma*cosh(i,j))+(0.5*a*senh(i,j))));
    end
end

%----Otra forma de Hallar A(t,T)
% q=sqrt((a^2)+(2*Var_in))
% 
% for i=1:N+1
%     for j=1:S
%         A1(i,j)=(2*a*b/Var_in)*log(2*q*exp((a+q)*((T-t(i,1))/2))...
%             /((a+q)*(exp(q*(T-t(i,1)))-1)+(2*q)));
%             
%     end
% end


%--------------Precio del Bono---------------------------------------------
for i=1:N+1
    for j=1:S
        B(i,j)=exp(A(i,j)-Rate(i,j)*D(i,j));
    end
end

%-------------Curva Rendimiento Vasicek Model------------------------------

for i=1:N+1
    for j=1:S
        R(i,j)=(-log(B(i,j)))/(T-t(i,1));
    end
end

%------------- Tasas Forward-----------------------------------------------

% for i=1:N+1
%     for j=1:S
%         F(i,j)= (b-(b-Rate(i,j)))*(exp(-a*(T-t(i,1))))...
%             -((Var/2)*D2(i,j));
%     end
% end

%----- Tasas promedio en T-------------------------------------------------

% TasaPromT=mean(Rate(N+1,:));

%----saca el promedio por simulacion

for i=1:N+1
        prom2(i,1)=mean(Rate(i,:));
end


%--- de prueba-------------------------------------------------------------
for i=1:N+1
        prom1(i,1)=mean(RR(i,:));
end

%--------------------------------------------------------------------------

% surf(Rate)
% title('CIR Model Short-Term')
% xlabel('Simulations')
% ylabel('Days')
% zlabel('Rate %')


%----Plotters--------------------------------------------------------------
subplot(3,3,2)
plot(0:dt:T,Rate)
title('r(t,T)-CIR- Moving')
subplot(3,3,1)
plot(0:dt:T,RR)
title('r(t,T)-CIR- Libre')
subplot(3,3,3)
plot(0:dt:T,D(:,1),'b')
title('D(T,T)=0')
subplot(3,3,4)
plot(0:dt:T,A(:,1),'r')
title('A(T,T)=0')
subplot(3,3,5)
plot(0:dt:T,B)
title('B(t,T)-CIR')
% subplot(3,3,5)
% plot(0:dt:T,R)
% title('R(t,T)-CIR')
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
title('Comparación Gráficas')



% legend('Tasa Real','Promedio')
% 
% subplot(3,3,10)
% boxplot(Rate), hold on
% plot(rt1)
% hold off
% title('CIR Model')
% xlabel('Simulations')

% plot(0:dt:T,rt1(:,1),'b',0:dt:T,prom1,'r',0:dt:T,prom2,'black')
% title('Comparación Gráficas - Modelo CIR')
% legend('Tasa Real','Simulación 1','Simulación 2')
% xlabel('Days')
% ylabel('Rates')