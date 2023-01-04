%% 1� Laborat�rio de Modela��o e Simula��o
% Alice Rosa, n�90007 
%
% Beatriz Pereira, n�90029 
%
% Grupo 16, Turno 3� feira �s 9h00 

 %% 1. Simula��o do movimento livre de uma viatura
 
%Cria��o do seguinte diagrama de blocos com o SIMULINK:

open('carro');

%%

%Defini��o dos par�metros de simula��o

m=30;
beta=5;
vo_set=[-3 3];
tau_set=[3 6 10];
yo=5;
tt=50;

close all;
f1=figure;
f2=figure;

for counter=1:length(vo_set);
    vo=vo_set(counter);
   
    for tau=[3 6 10]
    sim('carro',tt);
    
    figure(f1);
    plot(t,v);
    hold on;
    
    figure(f2);
    plot(t,y);
    hold on;
    end
end

figure(f1);
xlabel('t(s)','fontsize',12);
ylabel('v(m/s)','fontsize',12);
title('Velocidade','fontsize',12);
grid on;
set(gca,'fontsize',12);
lg=legend(sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(1)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(2)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(3)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(1)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(2)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(3)));
set(lg,'fontsize',9);   
hold on;

figure(f2);
xlabel('t(s)','fontsize',12);
ylabel('h(m)','fontsize',12);
title('Posi��o','fontsize',12);
grid on;
set(gca,'fontsize',12);
lg=legend(sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(1)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(2)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(1), tau_set(3)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(1)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(2)),...
       sprintf('V_o=%d m/s tau=%d s',vo_set(2), tau_set(3)));
set(lg,'fontsize',9);   

%%
% Tendo em conta que n�o est� a ser aplicada nenhuma for�a exterior sobre o
% carro, a �nica for�a que temos a atuar � a for�a de atrito, desta forma,
% � de esperar que passado algum tempo a velocidade do carro tenda para 0,
% que � exatamente o que se verifica nos gr�ficos obtidos. 
% Tamb�m podemos verificar que quanto menor a constante de tempo $\tau$,
% mais depressa o carro p�ra, o que seria de esperar tendo em conta a
% equa��o obtida na al�nea 1.3, pois a exponencial est� elevada a
% $-\frac{1}{\tau}$.
%
% Para os gr�ficos obtidos relativamente � posi��o do carro, esperava-se
% que passado algum tempo este estabilizasse numa posi��o que se pode
% calcular a partir da express�o obtida em 1.4: 
% $y_{final}=y_o+\frac{V_o*m}{\beta}$ 
%
% Tamb�m se verificou que para uma velocidade inicial negativa (carro em marcha atr�s)
% a posi��o em que o carro p�ra tamb�m � negativa tendo em conta o nosso
% referencial.

%% 2.Modelo Predador-Presa

%%% 2.2 Diferentes respostas para diferentes $\delta$ 

% Cria��o do seguinte diagrama de blocos com o SIMULINK:
open('predador_presa');
%%
clear;
close all;

delta1=1;
delta2=-1;
alpha1=1;
alpha2=1;
N1_inicial=3;
N2_inicial=3;
tt=30;

sim('predador_presa',tt);

figure(1);
plot(t,N1,t,N2,'LineWidth',1.2);
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
title('Solu��o para: \delta_1=1, \delta_2=-1','fontsize',12);
legend('N_1-Presa','N_2-Predador');
grid on;
hold on;

tt=3;
delta1=-1;
delta2=1;

sim('predador_presa',tt);
figure(2);
plot(t,N1,t,N2,'LineWidth',1.2);
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
title('Solu��o para: \delta_1=-1, \delta_2=1','fontsize',12);
legend('N_1-Presa','N_2-Predador');
grid on;
hold on;

tt=5;
delta1=-1;
delta2=-1;
N1_inicial=3;

sim('predador_presa',tt);
figure(3);
plot(t,N1,t,N2,'LineWidth',1.2);
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
title('Solu��o para: \delta_1=-1, \delta_2=-1','fontsize',12);
legend('N_1-Presa','N_2-Predador');
grid on;
hold on;


delta1=-1;
delta2=-4;

sim('predador_presa',tt);
figure(4);
plot(t,N1,t,N2,'LineWidth',1.2);
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
title('Solu��o para: \delta_1=-1, \delta_2=-4','fontsize',12);
legend('N_1-Presa','N_2-Predador');
grid on;
hold on;

%%
% Para o primeiro gr�fico, temos a solu��o oscilat�ria que como previsto
% ocorre quando $\delta_1$>0 e $\delta_2$<0. O segundo gr�fico � a
% solu��o em que o predador cresce indefenidamente e a presa se extingue,
% com $\delta_1<0$ e $\delta_2$>0. O terceiro gr�fico � a solu��o em que a
% presa e o predador se extiguem passado algum tempo, com $\delta_1$<0 e
% $\delta_2$<0. O �ltimo gr�fico obtido tem a mesma solu��o que o terceiro,
% no entanto, enquanto que no gr�fico anterior os predadores atingem um pico e
% depois,passado algum tempo, � que se extiguem, neste gr�fico diminui-se 
% o $\delta_2$ de forma a que as esp�cies se extiguissem sens�velmente ao
% mesmo tempo.

%% 2.3 Ponto de equil�brio, Espa�o de fase ($N_1$, $N_2$)

%%
clear;
close all;

N1_inicial_set=[3 2];
N2_inicial_set=[2 1];
delta1=1;
delta2=-1;
alpha1=1;
alpha2=1;
tt=30;

%Espa�o de fase

for counter=1:length(N1_inicial_set);
    N1_inicial=N1_inicial_set(counter);
    N2_inicial=N2_inicial_set(counter);
    
    sim('predador_presa',tt);
    
    figure(1);
    plot(N1,N2);
    hold on;
    
end

figure(1);
xlabel('Abund�ncia da Presa','fontsize',12);
ylabel('Abund�ncia do Predador','fontsize',12);
set(gca,'fontsize',12);
legend(sprintf('N1(0),N2(0)= %d,%d',N1_inicial_set(1),N2_inicial_set(1)),...
       sprintf('N1(0),N2(0)= %d,%d',N1_inicial_set(2),N2_inicial_set(2)));
title('Modelo predador-presa');
grid on;

%Ponto de equil�brio

N1_inicial=delta1/alpha1;
N2_inicial=-delta2/alpha2;
tt=10;

sim('predador_presa',tt);

figure(2);
p=plot(t,N1,t,N2,'--');
p(1).LineWidth = 1.5;
p(2).LineWidth = 2;
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
title('Ponto de equil�brio','fontsize',12);
legend('N1-Presa','N2-Predador');
grid on;

% Diferentes evolu��es temporais

N1_inicial_set=[1 2];
N2_inicial_set=[2 1];
tt=30;

for counter=1:length(N1_inicial_set);
    N1_inicial=N1_inicial_set(counter);
    N2_inicial=N2_inicial_set(counter);
    
    sim('predador_presa',tt);
    
    figure(3);
    p=plot(t,N1,t,N2);
    hold on;
end

figure(3);
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
set(gca,'fontsize',12);
lg=legend(sprintf('N_1-CI(%.0f,%.0f)',N1_inicial_set(1),N2_inicial_set(1)),...
       sprintf('N_2-CI(%.0f,%.0f)',N1_inicial_set(1),N2_inicial_set(1)),...
       sprintf('N_1-CI(%.1f,%.0f)',N1_inicial_set(2),N2_inicial_set(2)),...
       sprintf('N_2-CI(%.1f,%.0f)',N1_inicial_set(2),N2_inicial_set(2)));
set(lg,'fontsize',11);
title('Modelo predador-presa');
grid on;


%%
% No primeiro gr�fico obteve-se o espa�o de fase ($N_1$, $N_2$) para a solu��o
% oscilat�ria com diferentes valores de condi��es iniciais. 
% O segundo gr�fico � a confirma��o do ponto de equil�brio, ou seja,
% para determinadas condi��es iniciais (calculadas na al�nea 2.2) o sistema
% mant�m-se constante no tempo, que � o que se verifica. 
% O �ltimo gr�fico mostra como a mudan�a das condi��es iniciais alteram a
% evolu��o do sistema no tempo. A situa��o em que $N_1(0)$ < $N_2(0)$
% est� em avan�o e quadratura face � situa��o em que $N_1(0)$ > $N_2(0)$.

%% 2.4 Otimiza��o da curva

%%% 2.4 a)
 
%%
clear;
close all;

N1_inicial=4;
N2_inicial=1.6;
delta1=3.1;
delta2=-1.5;
alpha1=1.4;
alpha2=0.7;
tt=20;
load('presas.mat');

sim('predador_presa',tt);
plot(tr,yr,t,N2);
legend('N_1-presa','N_2-predador');
xlabel('Tempo','fontsize',12); ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
grid on;

%%
% O objectivo nesta quest�o era arranjar um equ�librio entre o $\alpha_2$ e 
% $N_2(0)$ de forma a quando o n�mero de presas diminuisse passado
% algum tempo o n�mero de predadores tamb�m come�ava a diminuir e
% vice-versa, mantendo sempre a mesma desfasagem entre as curvas da presa
% e predador. Desta forma, variou-se os par�metros at� se chegar a uma
% solu��o aproximada da que se queria sendo esta ($\alpha_2$, $N_2(0)$)=(0.7,1.6).

%% 2.4 b)

%%
clear;
close all;

N0_set=1.6:0.002:1.62;
alpha_set=0.7:0.002:0.72;
[X,Y]=meshgrid(N0_set,alpha_set);

for i=1:length(N0_set);
    V(1)=N0_set(i);
    
    for j=1:length(alpha_set);
        V(2)=alpha_set(j);
        dif(i,j)=MVAD(V);
    end
   % waitbar(i/10);
end

surf(X,Y,dif);
set(gca,'zscale','log');
set(gca,'zscale','log','fontsize',12);
xlabel('N_2(0)','fontsize',12); ylabel('\alpha_2','fontsize',12'); 
zlabel('Erro','fontsize',12); colorbar;

%%
% Fun��o utilizada para calcular o Erro

%%
type('MVAD.m');

%% 
% Nesta quest�o, o objectivo era encontrar os par�metros $\alpha_2$ e 
% $N_2(0)$ de forma a obtermos uma curva o mais pr�xima poss�vel da curva
% que nos foi dada. Para isto criou-se a fun��o erro que nos devolve o
% m�ximo valor absoluto das diferen�as entre os pontos fornecidos e os
% pontos simulados para determinados par�metros. No final, a partir do
% gr�fico obtido conclui-se que, dos par�metros testados, os que melhor se
% ajustam � curva fornecida s�o ($\alpha_2$, $N_2(0)$)=(0.704,1.6120).

%% 2.4 c)

clear;

N0_2=[1.6, 2];
alpha_2=[0.7, 3];

 for i=1:length(N0_2)
    xo=[N0_2(i),alpha_2(i)];
    fun=@MVAD;
    [point,erro]=fminsearch(fun,xo);
    display(erro);
    display(point);
 end

%%
% Nesta al�nea, para descobrir os par�metros utilizou-se a fun��o 
% fminsearch em que a partir da fun��o erro desenvolvida e de par�metros 
% iniciais esta fun��o vai fazendo chamadas � fun��o erro aproximando-se cada 
% vez mais do ponto desejado. Desta forma, o nosso c�digo fica muito mais 
% eficiente e os par�metros finais s�o mais exatos, relativamente � procura
% feita anteriormente, com ($\alpha_2$, $N_2(0)$)=(0.7047,1.6144). 
%
% No entanto, realizou-se ainda outro exemplo com outro ponto inicial 
% afastado do resultado ideal, onde se obteve um resultado completamente
% diferente. Podemos assim concluir que � melhor utilizar esta fun��o
% quando j� temos uma ideia dos par�metros a que corresponde o erro m�nimo,
% sen�o a fun��o encontra outro qualquer m�nimo local que n�o a resposta
% ideial.

%% 2.4 d)

clear;
close all;

N1_inicial=4;
N2_inicial=1.6144;
delta1=3.1;
delta2=-1.5;
alpha1=1.4;
alpha2=0.7047;
tt=20;
load('presas.mat');

sim('predador_presa',tt);
plot(tr,yr,'o',t,N1,'LineWidth',1.2);
legend('N1-Fornecido','N1-Simulado');
xlabel('Tempo','fontsize',12);
ylabel('Abund�ncia','fontsize',12);
set(gca,'fontsize',12);
legend('N_1-Pontos fornecidos','N_1-Curva simulada');
grid on;
ax = gca;
ax.GridAlpha = 0.6;

%%
% Esta al�nea serve como verifica��o dos resultados obtidos. Como podemos
% verificar utilizando os par�metros calculados na al�nea (c) a curva
% simulada ajusta-se muito bem aos pontos fornecidos.
 
%% 3. Sistema Ca�tico 
%%% 3.1. Realiza��o do modelo do sistema em Simulink e a passagem dos dados para o Matlab
%%
open('pendulo_duplo');
%%
clear;
close all;
 m=1;
 l=0.5;
 g=9.8;
 teta_1_inicial=0.05;
 teta_2_inicial=0.05;
 p_1_inicial=0;
 p_2_inicial=0;
 tt=3;
 
 sim('pendulo_duplo', tt);
 
 %posi��o ao longo do tempo 
x=l*(sin(teta_2)+sin(teta_1));
y=-l*(cos(teta_2)+cos(teta_1));

figure(1);
plot(x,y);
xlabel('x');
ylabel('y');
title('Sistema n�o ca�tico');

teta_1_inicial=15;
teta_2_inicial=15;
sim('pendulo_duplo', tt);

x=l*(sin(teta_2)+sin(teta_1));
y=-l*(cos(teta_2)+cos(teta_1));

figure(2);
plot(x,y);
xlabel('x');
ylabel('y');
title('Sistema ca�tico');

%%
% A partir dos gr�ficos obtidos verificou-se a implementa��o correta do
% sistema din�mico. Uma vez que para $\theta_1(0)$ e $\theta_2(0)$ pequenos
% obteve-se uma curva em que os �ngulos n�o se afastam muito dos �ngulos
% iniciais, como era de esperar. Enquanto que no segundo gr�fico, se
% verifica que para �ngulos maiores o sistema come�a a ter um
% comportamento mais aleat�rio e irregular.

%% 3.2. Curva de Lissajous
%%

clear;
close all;
 m=1;
 l=0.5;
 g=9.8;
 teta1_set=[0.05 5 10];
 teta2_set=[0.05 5 15];
 p_1_inicial=0;
 p_2_inicial=0;
 tt=2.5;
 titulos={'Curva de Lissajous - regime de fraca amplitude','Regime de m�dia amplitude',...
          'Regime de alta amplitude'};
 
for i=1:length(teta1_set)

teta_1_inicial=teta1_set(i);
teta_2_inicial=teta2_set(i);
sim('pendulo_duplo', tt);

figure(i);
plot(teta_1,teta_2);
title(titulos(i),'fontsize',14);
xlabel('\theta_1','fontsize',16);
ylabel('\theta_2','fontsize',16);
set(gca,'fontsize',14);
grid on;
end

%%
% Para valores muito pequenos de $\theta_1$ e $\theta_2$, verificamos que o
% sistema descreve uma curva de Lissajous quase perfeita, e que basta um
% pequeno aumento dos �ngulos para a curva come�ar a perder a sua forma. No
% �ltimo gr�fico com valores de $\theta$ inicial superiores a 10� a curva
% fica completamente desformada. Podemos concluir que dependendo dos
% �ngulos iniciais que escolhemos oferecemos menos/mais "momento" �s barras
% que originam traject�rias definidas/aleat�rias. 
%
% O tempo de simula��o influencia a forma da curva, uma vez que verific�mos
% que para, aproximadamente, tt=2.5 s o modelo desenha uma curva de
% Lissajous, enquanto que se aumentarmos o tempo(por exemplo, tt=10 s) o modelo desenha
% v�rias curvas.
% 
%% 3.4 Tempo at� que uma das barras fa�a um looping
%%

clear;
close all;
 m=1;
 l=0.5;
 g=9.8;
 tt=250;
 x=-1:0.05:1;
 y=-1:0.05:1;
 dteta1=0;
 dteta2=(-30*pi)/180;
 
 for i=1:length(x);
    
    for j=1:length(y);
        d=sqrt(x(i)^2+y(j)^2);
      %Verifica��o se a dist�ncia � inferior a 2*l
      if d > 2*l
          tempo(i,j)=NaN;
          continue;
      else
        %Encontrar a posi��o do pend�lo atrav�s das coordenadas
        teta=posicao(x(i),y(j),l);
        
        teta_1_inicial=teta(1);
        teta_2_inicial=teta(2);
        p_1_inicial=(1/6)*m*(l^2)*(3*dteta2*cos(teta_1_inicial-teta_2_inicial));
        p_2_inicial=(1/6)*m*(l^2)*(2*dteta2);
        
        sim('pendulo_duplo',tt);
        %Encontrar o ind�ce em que o p�ndulo d� um loop (volta de 360�)
        k=find((teta_1<teta_1_inicial-(2*pi)) | (teta_1>teta_1_inicial+(2*pi)) | ...
               (teta_2<teta_2_inicial-(2*pi)) | (teta_2>teta_2_inicial+(2*pi)),1);
        
        if k~=0 %Se k=0, nenhuma das barras deu um loop
            tempo(i,j)=t(k);
        else
            tempo(i,j)=NaN;
        end
        
      end
   end
    %waitbar(i/(length(x)));
 end
 
figure(1);
pcolor(x,y,log10(tempo'));
xlabel('x','fontsize',12);
ylabel('y','fontsize',12);
set(gca,'fontsize',12);
title('Tempos de looping','fontsize',12);
colorbar;

%% 
% Para a realiza��o desta al�nea foi utilizada a fun��o 'posicao', que
% recebe as coordenadas da ponta do p�ndulo e retorna os �ngulos iniciais
% do mesmo, utilizando as express�es obtidas na al�nea 3.3.

%%
type('posicao.m');

%%
% A partir dos resultados obtidos verificou-se que para y mais baixos, maior
% � o tempo necess�rio para que uma das barras d� um looping, uma vez que
% quanto menor a coordenada y, menor � a energia pontencial na ponta do p�ndulo, 
% logo, � necess�rio uma velocidade angular inicial superior para levar a uma 
% das barras do p�ndulo a realizar um loop.

 %% 3.4 Verifica��o de resultados
 %%
  clear;
 close all;
 m=1;
 l=0.5;
 g=9.8;
 tt_set=[30 100 250];
 intervalos=[0 30 100 250];
 x=[-0.7,-0.45,0.1];
 y=[-0.4,-0.65,-0.55];
 dteta1=0;
 dteta2=(-30*pi)/180;
 titulos={'Ponto de looping t=[0 30]s','Ponto de looping t=[30 100]s',...
          'Ponto de looping t=[100 250]s'};
 
 for i=1:length(x);
       tt=tt_set(i);
       teta=posicao(x(i),y(i),l);
         
       teta_1_inicial=teta(1);
       teta_2_inicial=teta(2);
       p_1_inicial=(1/6)*m*(l^2)*(3*dteta2*cos(teta_1_inicial-teta_2_inicial));
       p_2_inicial=(1/6)*m*(l^2)*(2*dteta2);
        
       sim('pendulo_duplo',tt);
         
       k=find((teta_1<teta_1_inicial-(2*pi)) | (teta_1>teta_1_inicial+(2*pi)) | ...
              (teta_2<teta_2_inicial-(2*pi)) | (teta_2>teta_2_inicial+(2*pi)),1);
       
       figure(i);
       plot(t,teta_2);
       hold on;
       ponto=plot(t(k),teta_2(k),'r*');
       xlim([intervalos(i) intervalos(i+1)]);
       title(titulos(i),'fontsize',12);
       xlabel('Tempo [s]','fontsize',12); ylabel('\theta_2 [rad]','fontsize',12);
       set(gca,'fontsize',12);
       legend(ponto,sprintf('Ponto de looping=(%.2f s,%.2f rad)',t(k),teta_2(k)));
       grid on;
       
       figure(i+3);
       plot(t,teta_1);
       xlim([intervalos(i) intervalos(i+1)]);
       title(titulos(i),'fontsize',12);
       xlabel('Tempo [s]','fontsize',12); ylabel('\theta_1 [rad]','fontsize',12);
       set(gca,'fontsize',12);
       grid on;
 end
 
 %%
 % Tendo em considera��o que a primeira barra tem uma velocidade angular nula, � muito
 % pouco prov�vel que esta barra seja a primeira a dar o loop. Desta forma,
 % neste caso, verificou-se que apenas a segunda barra d� o loop.
 