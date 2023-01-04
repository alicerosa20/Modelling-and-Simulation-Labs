%% 4� Laborat�rio de Modela��o e Simula��o 2019/20
% Detec��o de hotspots WiFi
%
% Alice Rosa, N� 90007 
%
% Beatriz Pereira, N� 90029 
%
% Grupo 16, Turno 3� feira �s 9h00 
%% Quest�o 2.a)
%
% Nesta quest�o, decomp�s-se a matriz transposta de probabilidades de transi��o
% entre estados, P', em valores pr�prios e vetores pr�prios. De seguida, obteve-se a 
% distribui��o de equil�brio da cadeia de Markov
% por normaliza��o do vetor pr�prio associado ao valor pr�prio mais pr�ximo
% de 1.
%%
clear;
close all;

load('MarkovChain.mat');

%valores e vetores pr�prios da matriz P'
[vect, val]= eig(P');

%Determinar o valor pr�prio mais pr�ximo de 1
min= 900;

  for i=1:length(val) 
    if(abs(val(i,i)-1))<min
        min=abs(val(i,i)-1);
        ind=i;
    end
  end
% Associar o vetor pr�prio a esse �ndice
Vp= vect(:,ind);

% normaliza��o do vetor pr�prio
soma_Vp= sum(Vp);
Vp_normalizado = Vp/soma_Vp;

%soma das probabilidades tem de ser igual a 1
Vp_normalizado_total = sum(Vp_normalizado);

%distribui��o de equilibro da cadeia de Markov
bar(Vp_normalizado);

xlabel('Estados');
xlim([0 21]);
ylabel('Probabilidade de estado');
title('Distribui��o de equil�brio da cadeia de Markov');
set(gca,'Fontsize',12);
%%
% Por an�lise do gr�fico da distribui��o de equil�brio da cadeira de
% Markov, conclui-se que os estados mais prov�veis s�o o 7 e o 19
% (P=0.09649), sendo que estes correspondem, tamb�m, aos n�s do grafo com
% maior n�mero de liga��es e os estados menos prov�veis s�o o 8 e 17 (P=0.01072). 

%% Quest�o 2.b)
%
% Tendo em conta a distribui��o de equil�brio da cadeia de Markov da al�nea
% anterior, resolveu-se uma vers�o ponderada do problema de m�nimos 
% quadr�ticos, para um M suficientemente longo de medi��es arbitr�rio, 
% de modo a determinar a posi��o da fonte estimada.
%%
Po=100; %Pot�ncia da fonte

variancia=10^-2;

M=1000; %N�mero de medi��es

a=[nodePos(:,2),nodePos(:,3)]';

D=squareform(pdist([sourcePos' zeros(size(sourcePos')) a]'));
d=D(1,3:end); %Fonte-�ncora dist�ncia

% Posi��o estimada da fonte
x=sourcePos_estimativa(M,Po,variancia,Vp_normalizado,a,d);

% Plot 
figure(2)
plot(nodePos(:,2),nodePos(:,3),'o'); hold on
plot(sourcePos(:,1),sourcePos(:,2),'x'); plot(x(:,1),x(:,2),'s'); hold on
axis(100*[0 1 0 1]); axis('square')
legend( {'�ncoras','Fonte', 'Posi��o estimada da fonte'},'Location','southwest' );
title('Posi��o da fonte e das �ncoras');

% Fun��o que estima a posi��o da fonte
type ('sourcePos_estimativa');
%%
% Conclui-se que, para um n�mero elevado de medi��es, o m�todo dos m�nimos
% quadrados fornece uma boa estimativa para a posi��o da fonte, sendo que o
% pequeno desvio apresentado deve-se � presen�a de ru�do.

%% Quest�o 2.c)
%
% Fez-se a simula��o da evolu��o das probabilidades dos diversos estados da cadeia de 
% Markov, ao longo do tempo, para diferentes condi��es iniciais de $\pi$(0).
%
% A evolu��o das probabilidades dos estados de uma cadeia de Markov,
% entre dois instantes de tempo consecutivos, � dada por
% $\pi'(t+1)=\pi'(t)P$.
%%

clearvars -except Vp_normalizado;
close all;

load('MarkovChain.mat');

% --- 1� Condi��o Inicial ---

 t=1:200; 
 
 pi_0=zeros(1,20)';
 
 % Inicializa��o no estado 7
 pi_0(7)=1;
 pi_t1=pi_0;
 
 %soma de todas as entradas do vetor pi
 soma_pi_1(1)=sum(pi_t1);
 
 % evolu��o das probabilidades dos estados de uma cadeia de Markov
 for i=1:199
     pi_t1(:,i+1)=pi_0'*(P^(i));
     soma_pi_1(i+1)=sum(pi_t1(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
figure(1)
plot3(t,estados,pi_t1);
xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Condi��o inicial no centro do grafo');
set(gca,'Fontsize',12);
 
 
 figure(2) 
 plot(t,soma_pi_1);
 ylim([0.95 1.05]);
 xlabel('Tempo (t)');
 ylabel('Soma das probabilidades');
 title('Soma das entradas \pi(t)');
 set(gca,'Fontsize',12);
 
 % --- 2� Condi��o Inicial ---
 
 pi_0=zeros(1,20)';
 
 % Inicializa��o no estado 17
 pi_0(17)=1;
 pi_t2=pi_0;
 
 soma_pi_2(1)=sum(pi_t2); %soma das probabilidades
 
 % evolu��o das probabilidades dos estados de uma cadeia de Markov
 for i=1:199
     pi_t2(:,i+1)=pi_0'*(P^(i));
     soma_pi_2(i+1)=sum(pi_t2(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
figure(3)
plot3(t,estados,pi_t2);
xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Condi��o inicial num estado que ret�m o token');
set(gca,'Fontsize',12);
 
 
 figure(4) 
 plot(t,soma_pi_2);
 ylim([0.95 1.05]);
 xlabel('Tempo (t)');
 ylabel('Soma das probabilidades');
 title('Soma das entradas \pi(t)');
 set(gca,'Fontsize',12);
 
 % --- 3� Condi��o Inicial ---
 
 pi_0=zeros(1,20)';
 %distribui��o inicial uniforme 
 pi_0(:,1)=1/20;
 pi_t3=pi_0;
 
 soma_pi_3(1)=sum(pi_t3); %soma das probabilidades

 % evolu��o das probabilidades dos estados de uma cadeia de Markov
 for i=1:199
     pi_t3(:,i+1)=pi_0'*(P^(i));
     soma_pi_3(i+1)=sum(pi_t3(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
figure(5)
plot3(t,estados,pi_t3);
xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Distribui��o inicial uniforme');
set(gca,'Fontsize',12);
 
 figure(6) 
 plot(t,soma_pi_3);
 ylim([0.95 1.05]);
 xlabel('Tempo (t)');
 ylabel('Soma das probabilidades');
 title('Soma das entradas \pi(t)');
 set(gca,'Fontsize',12);
 
 % --- 4� Condi��o Inicial ---
 
 %distribui��o inicial em equilibrio
 pi_0= Vp_normalizado;
 pi_t4=pi_0;
 
 soma_pi_4(1)=sum(pi_t4); %soma das probabilidades
 
 for i=1:199
     pi_t4(:,i+1)=pi_0'*(P^(i));
     soma_pi_4(i+1)=sum(pi_t4(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
figure(7)
plot3(t,estados,pi_t4);
xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Distribui��o inicial em equil�brio');
set(gca,'Fontsize',12);
 
 
 figure(8) 
 plot(t,soma_pi_4);
 ylim([0.95 1.05]);
 xlabel('Tempo (t)');
 ylabel('Soma das probabilidades');
 title('Soma das entradas \pi(t)');
 set(gca,'Fontsize',12); 
 %%
 % Simulou-se a evolu��o das probabilidades dos diversos estados da cadeia de 
 % Markov, ao longo do tempo, para quatro condi��es iniciais diferentes de $\pi(0)$.
 % 
 % Em primeiro lugar, incializou-se o token na �ncora 7, pois esta 
 % � uma �ncora central do grafo com quatro liga��es todas com igual probabilidade (P=0.25).
 % Verifica-se que as probabilidades de todos os estados tendem para o equil�brio.
 %
 % De seguida, inicializou-se o token num dos subconjuntos em que este fica
 % retido durante um longo per�odo de tempo por ter menos probabilidade de sair
 % do que ficar, portanto escolheu-se a �ncora 17, que
 % corresponde, tamb�m, a um dos estados com menor probabilidade. Uma vez que
 % se inicializa o token dentro deste subconjunto e este s� t�m 20% 
 % probabilidade de sair, verifica-se, por observa��o do gr�fico, que o tempo de estabelecimento � 
 % superior ao caso anterior, isto �, as probabilidades dos estados demoram 
 % mais tempo a atingir o equil�brio. 
 %
 % Como terceiro caso, escolheu-se um vetor de probabilidades com
 % distribui��o inicial uniforme. Daqui, uma vez que o token tem a mesma probabilidade
 % de come�ar em qualquer uma das �ncoras, este corresponde ao caso m�dio dos
 % dois anteriores.
 %
 % Por �ltimo, escolheu-se um vetor de probabilidades com
 % distribui��o inicial em equilibrio, em que este � igual ao vetor de 
 % vetores pr�prios correspondente �
 % distribui��o de equil�bro da cadeia de Markov (decomposi��o em valores 
 % e vetores pr�prios). Uma vez que o sistema se mant�m 
 % nesse estado indefinidamente, comprova-se que a distribui��o do sistema 
 % tende sempre para a posi��o de equil�brio.
 
 %% Quest�o 2.d)
 clear;
 close all;
 
 MarkovChainDraw
 title('Grafo de comunica��es entre os agentes');
 %%
 % Atrav�s da an�lise do grafo de comunica��es entre agentes, das 
 % respetivas probabilidades de transi��o entre os v�rios estados e da
 % distribui��o de equil�brio, identificaram-se dois subconjuntos de
 % estados da cadeia onde o token pode ficar a circular durante relativamente 
 % longos intervalos de tempo.
 %
 % Subconjunto 1: 5, 6, 11, 15.
 %
 % Subconjunto 2: 8, 9, 10, 12, 17.
 %
 % O token tem tend�ncia a ficar preso nestes subconjuntos durante per�odos
 % mais longos, pois em ambos os casos s� h� uma forma de sa�da e esta
 % apresenta uma probabilidade inferiror �s liga��es de entrada. O
 % subconjunto 1, apresenta uma probabilidade de sair de 20%, apenas a
 % partir do n� 6, enquanto a probabilidade de se manter � de 80%. O mesmo
 % acontece no subconjunto 2, em o token apresenta 20% de probabilidade de
 % sair pelo n� 12 e 80% de ficar e se deslocar para o n� 8 ou 10.
 %
 % Tal como se verificou na al�nea anterior, se o token for inicializado
 % num destes subconjuntos ou entrar num deles, ir� tender para o equil�brio, 
 % mas o seu tempo de converg�ncia ser� maior.
 %% Quest�o 2.d) Ponto de vista da equipa
 % Para melhorar a efic�cia de circula��o global do token e evitar que
 % este fique preso nestes subconjuntos, procedeu-se � altera��o dos pesos dos n�s 1, 3,
 %  4, 6, 10, 12, 13, 15, 16 e 20 e adicionaram-se liga��es estrat�gicas
 %  entre os n�s 3 e 16, 15 e 20, 13 e 1.
 %%
close all 
clear

load('MarkovChain.mat');

% Altera��o dos pesos dos n�s na matriz P
P(6,1)=0.3; P(6,11)=0.45; P(6,15)=0.25; %N� 6
P(15,20)=0.25; P(15,5)=0.40; P(15,6)=0.35; %N� 15
P(1,6)=0.20; P(1,20)=0.25; P(1,7)=0.25; P(1,13)=0.30; %N� 1
P(3,12)=0.4; P(3,19)=0.3; P(3,16)=0.3; %N� 3
P(12,3)=0.30; P(12,10)=0.25; P(12,8)=0.45; %N� 12
P(4,19)=0.4; P(4,13)=0.15; P(4,2)=0.45; %N� 4
P(20,15)=0.35; P(20,1)=0.15; P(20,7)=0.15; P(20,14)=0.35; %N� 20
P(13,1)=0.25; P(13,2)=0.35; P(13,4)=0.15; P(13,19)=0.25; %N� 13
P(16,3)=0.25; P(16,7)=0.15; P(16,18)=0.6; %N� 16
P(10,17)=0.45; P(10,12)=0.35; P(10,9)=0.2; %N� 10

% repetir al�nea 2.a)

%valores e vetores pr�prios da matriz P'
[vect, val]= eig(P');

%Determinar o valor pr�prio mais pr�ximo de 1
min= 900;

  for i=1:length(val) 
    if(abs(val(i,i)-1))<min
        min=abs(val(i,i)-1);
        ind=i;
    end
  end

Vp= vect(:,ind);

% normaliza��o 
soma_Vp= sum(Vp);
Vp_normalizado = Vp/soma_Vp;

%soma das probabilidades tem de ser igual a 1
Vp_normalizado_total = sum(Vp_normalizado);

%distribui��o de equilibro da cadeia de Markov
bar(Vp_normalizado);

xlabel('Estados');
xlim([0 21]);
ylabel('Probabilidade de estado');
title('Distribui��o de equil�brio da cadeia de Markov melhorada');
set(gca,'Fontsize',12);

%repetir 2.b)

Po=100; %Pot�ncia da fonte

variancia=10^-2;

M=1000; %N�mero de medi��es

a=[nodePos(:,2),nodePos(:,3)]';

D=squareform(pdist([sourcePos' zeros(size(sourcePos')) a]'));
d=D(1,3:end); %Fonte-�ncora dist�ncia

% Estimativa da localiza��o da fonte
x=sourcePos_estimativa(M,Po,variancia,Vp_normalizado,a,d);

% Plot %
figure
plot(nodePos(:,2),nodePos(:,3),'o'); hold on
plot(sourcePos(:,1),sourcePos(:,2),'x'); plot(x(:,1),x(:,2),'s'); hold on
axis(100*[0 1 0 1]); axis('square')
legend( {'�ncoras','Fonte', 'Posi��o estimada da fonte'},'Location','southwest' );
title('Posi��o da fonte e das �ncoras com grafo melhorado');

% repetir 2.c)

t=1:200; 
 
 pi_0=zeros(1,20)';
 pi_0(7)=1;
 pi_t1=pi_0;
 soma_pi_1(1)=sum(pi_t1);
 
 for i=1:199
     pi_t1(:,i+1)=pi_0'*(P^(i));
     soma_pi_1(i+1)=sum(pi_t1(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
 figure
 plot3(t,estados,pi_t1);
 xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Grafo Melhorado');
set(gca,'Fontsize',12);
 

%% Quest�o 2.d) Ponto de vista de um elemento hostil
% Fizeram-se altera��es no peso das liga��es dos n�s 1 e 6 (jamming
% selectivo do canal de comunica��es), de modo dificultar a circula��o 
% do token e aumentar o seu isolamento nos subconjuntos que o ret�m durante
% per�odos de tempo mais longos.
%%
close all
clear

load('MarkovChain.mat');

% Altera��o na matriz P
P(1,6)=0.8; P(1,20)=0.1; P(1,7)=0.1;
P(6,1)=0.1; P(6,11)=0.45; P(6,15)=0.45;

% repetir al�nea 2.a)

%valores e vetores pr�prios da matriz P'
[vect, val]= eig(P');

%Determinar o valor pr�prio mais pr�ximo de 1
min= 900;

  for i=1:length(val) 
    if(abs(val(i,i)-1))<min
        min=abs(val(i,i)-1);
        ind=i;
    end
  end

Vp= vect(:,ind);

% normaliza��o 
soma_Vp= sum(Vp);
Vp_normalizado = Vp/soma_Vp;

%soma das probabilidades tem de ser igual a 1
Vp_normalizado_total = sum(Vp_normalizado);

%distribui��o de equilibro da cadeia de Markov
bar(Vp_normalizado);

xlabel('Estados');
xlim([0 21]);
ylabel('Probabilidade de estado');
title('Distribui��o de equil�brio da cadeia de Markov hostil');
set(gca,'Fontsize',12);

%repetir 2.b)

Po=100; %Pot�ncia da fonte

variancia=10^-2;

M=1000; %N�mero de medi��es

a=[nodePos(:,2),nodePos(:,3)]';

D=squareform(pdist([sourcePos' zeros(size(sourcePos')) a]'));
d=D(1,3:end); %Fonte-�ncora dist�ncia

x=sourcePos_estimativa(M,Po,variancia,Vp_normalizado,a,d);

% Plot %
figure
plot(nodePos(:,2),nodePos(:,3),'o'); hold on
plot(sourcePos(:,1),sourcePos(:,2),'x'); plot(x(:,1),x(:,2),'s'); hold on
axis(100*[0 1 0 1]); axis('square')
legend( {'�ncoras','Fonte', 'Posi��o estimada da fonte'},'Location','southwest' );
title('Posi��o da fonte e das �ncoras com elemento hostil');

% repetir 2.c)

t=1:800; 
 
 pi_0=zeros(1,20)';
 pi_0(7)=1;
 pi_t1=pi_0;
 soma_pi_1(1)=sum(pi_t1);
 
 for i=1:799
     pi_t1(:,i+1)=pi_0'*(P^(i));
     soma_pi_1(i+1)=sum(pi_t1(:,i+1));
 end
 
 estados=repmat(1:20,length(t),1);
 
 figure
 plot3(t,estados,pi_t1);
  xlabel('Tempo (t)');
ylabel('Estados');
zlabel('Probabilidade de cada estado');
title('Grafo hostil');
set(gca,'Fontsize',12);

%%
% Conclu�-se, destas altera��es, que a fluidez de circula��o do token tem  
% influ�ncia directa na precis�o de localiza��o da fonte, uma vez que ao
% adicionar o elemento hostil, o tempo que a distribui��o de probabilidades 
% demora a tender para a distribui��o de equil�brio � muito mais elevado em 
% compara��o com o do o grafo melhorado e a localiza��o da fonte � muito
% menos precisa do que quando ocorre uma circula��o eficaz do token.

 %% Quest�o 3.a)
 % Nesta quest�o, introduz-se o m�todo de Monte Carlo, onde se simula o
 % avan�o do token que parte de um estado inicial e avan�a aleatoriamente
 % pelo diagrama de estados durante um certo n�mero de instantes de tempo,
 % anotando-se o historial. Este procedimento, designado por run de Monte
 % Carlo � repetido um n�mero suficiente de vezes de forma a assegurar que
 % os resultados t�m significado estat�stico. 
 %
 % Foi desenvolvida a fun��o 'simulacao_MonteCarlo' de forma a simular-se
 % para um certo n�mero de runs (n_runs) e de passos (n_passos) o avan�o do
 % token pelo diagrama de estados e guardar essa informa��o para futura
 % an�lise.
 %%
 type('simulacao_MonteCarlo.m')
 %%
 
 clear 
 close all
 
 load('MarkovChain.mat')
 
 n_runs=5000; % n�mero total de runs
 n_passos=300; % n�mero total de passos por run
 
 estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,20,'a',0);
 %n�mero total de passos por estado 
 estados_MC = cumsum(estados_MC);
 
 estados_MC_norm=zeros(n_passos, 20);
 %normaliza��o da matriz de estados por passo
 for j=1:n_passos
     estados_MC_norm(j,:)=estados_MC(j,:)./sum(estados_MC(j,:));
 end
 %Distribui��o de equ�librio das probabilidades encontra-se na �ltima
 %linha 
 dist_estados=estados_MC_norm(n_passos,:);
 
 %repeti��o da al�nea 2.a)
 
 [vect, val]= eig(P');
 
 min= 900;
 
 for i=1:length(val)
     if(abs(val(i,i)-1))<min
         min=abs(val(i,i)-1);
         ind=i;
     end
 end
 
 Vp= vect(:,ind);
 
 soma_Vp= sum(Vp);
 Vp_normalizado = Vp/soma_Vp;
 
 Vp_normalizado_total = sum(Vp_normalizado);

%Compara��o das distribui��es de equil�brio 
 figure;
 bar([Vp_normalizado dist_estados']);
 xlabel('Estados'); xlim([0 21]);
 ylabel('Probabilidade de estado'); title('Distribui��es de equil�brio');
 legend('Original','Monte Carlo');
 
 %Evolu��o no tempo da distribui��o
 t=1:n_passos;
 estados=repmat(1:20,length(t),1);
 
 figure;
 plot3(t,estados,estados_MC_norm);
 xlabel('t [s]'); ylabel('�ncoras'); zlabel('Probabilidade');
 title('Evolu��o no tempo da distribui��o de probabilidade');
 
 %-----Compara��o com o estado mais frequente----- 
 
 estado_inicial=7;
 %Obten��o da matriz passo/estado
 estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,20,'c',estado_inicial);
 
 estados_MC = cumsum(estados_MC);
 
 estados_MC_norm=zeros(n_passos, 20);
 
 for j=1:n_passos
     estados_MC_norm(j,:)=estados_MC(j,:)./sum(estados_MC(j,:));
 end
 
 dist_estados=estados_MC_norm(n_passos,:);

%Compara��o das distribui��es de equil�brio 
 figure;
 bar([Vp_normalizado dist_estados']);
 xlabel('Estados'); xlim([0 21]);
 ylabel('Probabilidade de estado'); 
 title('Distribui��es de equil�brio - Token inicial mais frequente');
 legend('Original','Monte Carlo');
 
 %Evolu��o no tempo da distribui��o
 figure;
 plot3(t,estados,estados_MC_norm);
 xlabel('t [s]'); ylabel('�ncoras'); zlabel('Probabilidade');
 title('Evolu��o da distribui��o no tempo - Token inicial mais frequente');
 
  %-----Compara��o com o estado menos frequente----- 
 
 estado_inicial=17;
 
 n_runs=5000;
 n_passos=300;
 
 estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,20,'c',estado_inicial);
 
 estados_MC = cumsum(estados_MC);
 
 estados_MC_norm=zeros(n_passos, 20);
 
 for j=1:n_passos
     estados_MC_norm(j,:)=estados_MC(j,:)./sum(estados_MC(j,:));
 end
 
 dist_estados=estados_MC_norm(n_passos,:);

%Compara��o das distribui��es de equil�brio 
 figure;
 bar([Vp_normalizado dist_estados']);
 xlabel('Estados'); xlim([0 21]);
 ylabel('Probabilidade de estado'); 
 title('Distribui��es de equil�brio - Token inicial menos frequente');
 legend('Original','Monte Carlo');
 
 %Evolu��o no tempo da distribui��o
 figure;
 plot3(t,estados,estados_MC_norm);
 xlabel('t [s]'); ylabel('�ncoras'); zlabel('Probabilidade');
 title('Evolu��o da distribui��o no tempo - Token inicial menos frequente');
 
 %%
 % A partir da primeira figura, podemos confirmar que as distribui��es de
 % equil�brio obtidas pelo m�todo de Monte Carlo s�o pr�ximas das
 % determinadas na sec��o 2.
 %
 % De forma semelhante � al�nea 2.c), a partir do
 % m�todo de Monte Carlo, simulou-se a evolu��o das
 % probabilidades dos diversos estados da cadeira de Markov ao longo do
 % tempo para diferentes condi��es iniciais. 
 %
 % Pode-se observar que, quando o token come�a num estado central,
 % como o estado 7, este tem mais facilidade em rapidamente circular pela
 % rede e atingir a probabilidade de equil�brio, ou seja, tem um elevado 
 % ritmo de converg�ncia. O mesmo n�o se verifica quando o token inicia
 % numa zona que o ret�m, por exemplo, o estado 17. Neste caso, o ritmo de
 % converg�ncia � bastante inferior e as probabilidades de equil�brio
 % atingidas s�o diferentes das originais.
 
 %%
 % 3.a) Ponto de vista da equipa/elemento hostil
 %%
 clear 
 
 load('MarkovChain.mat')
 
 P(6,1)=0.3; P(6,11)=0.45; P(6,15)=0.25; %N� 6
 P(15,20)=0.25; P(15,5)=0.40; P(15,6)=0.35; %N� 15
 P(1,6)=0.20; P(1,20)=0.25; P(1,7)=0.25; P(1,13)=0.30; %N� 1
 P(3,12)=0.4; P(3,19)=0.3; P(3,16)=0.3; %N� 3
 P(12,3)=0.30; P(12,10)=0.25; P(12,8)=0.45; %N� 12
 P(4,19)=0.4; P(4,13)=0.15; P(4,2)=0.45; %N� 4
 P(20,15)=0.35; P(20,1)=0.15; P(20,7)=0.15; P(20,14)=0.35; %N� 20
 P(13,1)=0.25; P(13,2)=0.35; P(13,4)=0.15; P(13,19)=0.25; %N� 13
 P(16,3)=0.25; P(16,7)=0.15; P(16,18)=0.6; %N� 16
 P(10,17)=0.45; P(10,12)=0.35; P(10,9)=0.2; %N� 10

 n_runs=5000;
 n_passos=300;
 
 estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,20,'a',0);
 
 estados_MC = cumsum(estados_MC);
 
 estados_MC_norm=zeros(n_passos, 20);
 
 for j=1:n_passos
     estados_MC_norm(j,:)=estados_MC(j,:)./sum(estados_MC(j,:));
 end
 
 dist_estados=estados_MC_norm(n_passos,:);
 
 %repeti��o da al�nea 2.a)
 
 [vect, val]= eig(P');
 
min= 900;

  for i=1:length(val) 
    if(abs(val(i,i)-1))<min
        min=abs(val(i,i)-1);
        ind=i;
    end
  end

Vp= vect(:,ind);

soma_Vp= sum(Vp);
Vp_normalizado = Vp/soma_Vp;
Vp_normalizado_total = sum(Vp_normalizado);

%Compara��o das distribui��es de equil�brio 
 figure;
 bar([Vp_normalizado dist_estados']);
 xlabel('Estados'); xlim([0 21]);
 ylabel('Probabilidade de estado'); 
 title('Distribui��es de equil�brio - Grafo melhorado');
 legend('Original','Monte Carlo');
 
 %Evolu��o no tempo da distribui��o
 t=1:n_passos;
 estados=repmat(1:20,length(t),1);
 
 figure;
 plot3(t,estados,estados_MC_norm);
 xlabel('t [s]'); ylabel('�ncoras'); zlabel('Probabilidade');
 title('Evolu��o da distribui��o no tempo - Grafo melhorado');
 
 % 3.a) Ponto de vista de um elemento hostil 
 
 clear P
 load('MarkovChain.mat')
 
 P(1,6)=0.8; P(1,20)=0.1; P(1,7)=0.1;
 P(6,1)=0.1; P(6,11)=0.45; P(6,15)=0.45;
 
 estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,20,'a',0);
 
 estados_MC = cumsum(estados_MC);
 
 estados_MC_norm=zeros(n_passos, 20);
 
 for j=1:n_passos
     estados_MC_norm(j,:)=estados_MC(j,:)./sum(estados_MC(j,:));
 end
 
 dist_estados=estados_MC_norm(n_passos,:);
 
 %al�nea 2.a)
 
 [vect, val]= eig(P');

 min= 900;
 
 for i=1:length(val)
     if(abs(val(i,i)-1))<min
         min=abs(val(i,i)-1);
         ind=i;
     end
 end

 Vp= vect(:,ind);
 
 soma_Vp= sum(Vp);
 Vp_normalizado = Vp/soma_Vp;
 Vp_normalizado_total = sum(Vp_normalizado);

%Compara��o das distribui��es de equil�brio 
 figure;
 bar([Vp_normalizado dist_estados']);
 xlabel('Estados'); xlim([0 21]);
 ylabel('Probabilidade de estado'); 
 title('Distribui��es de equil�brio - Grafo hostil');
 legend('Original','Monte Carlo');
 
 %Evolu��o no tempo da distribui��o
 figure;
 plot3(t,estados,estados_MC_norm);
 xlabel('t [s]'); ylabel('�ncoras'); zlabel('Probabilidade');
 title('Evolu��o da distribui��o no tempo - Grafo hostil');
 
 %%
 % A partir do m�todo de Monte Carlo (MC) tamb�m se simulou as variantes do grafo da
 % al�nea 2.d). Para a cadeia de Markov melhorada a
 % distribui��o de equil�brio obtida � ainda mais pr�xima da original do que a
 % determinada anteriormente, tal pode ser explicado pelo facto de nesta
 % o vetor de probabilidade limite ser pr�ximo de uma distribui�ao
 % uniforme.
 %
 % Relativamente aos ritmos de converg�ncia, conclui-se que estes s�o mais
 % lentos para o m�todo de MC. Para a cadeia de Markov com pior circula��o,
 % � bastante mais lento pois nem chega a atingir a distribui��o de
 % equil�brio dentro da janela de tempo.
 %
 % Desta forma, pode-se concluir que o m�todo de MC � bastante preciso para
 % uma quantidade elevada de runs.
 %% Quest�o 3.b)
 % Nesta quest�o, estimou-se o erro da posi��o da fonte ao longo do tempo a
 % partir da fun��o 'erro_MonteCarlo.m'. Esta constr�i as matrizes A e
 % b, resolve o problema de m�nimos quadr�ticos a partir do algoritmo RLS e 
 % calcula o erro atrav�s da posi��o estimada da fonte em cada passo.
 %%
 type('erro_MonteCarlo.m')
 %%
 clear
 close all
 
 load('MarkovChain.mat')
 n_passos=200;
 n_runs=800;

 Po=100; %Pot�ncia da fonte
 variancia=10^-2;
 
 estado_inicial=7;
 
 a=[nodePos(:,2),nodePos(:,3)]';
 D=squareform(pdist([sourcePos' zeros(size(sourcePos')) a]'));
 d=D(1,3:end); %Fonte-�ncora dist�ncia
 
 erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,d,P,estado_inicial,'c','a',0);

 erro_med=erro./n_runs;
 t=1:n_passos;
 
 figure(1)
 plot(t, erro_med);
 hold on
 
 % Zona 1 longe da fonte
 
 estado_inicial=15;
 erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,d,P,estado_inicial,'c','a',0);

 erro_med=erro./n_runs;
 t=1:n_passos;
 
 figure(1)
 plot(t, erro_med);
 hold on
 
 % Zona 2 perto da fonte
 
  estado_inicial=10;
 erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,d,P,estado_inicial,'c','a',0);

 erro_med=erro./n_runs;
 t=1:n_passos;
 
 figure(1)
 plot(t, erro_med);
 xlabel('t [s]'); ylabel('Erro');
 title('Erro de estimativa de posi��o da fonte');
 legend('Posi��o Central','Zona 1 - Longe da Fonte', 'Zona 2 - Perto da fonte');

 % Vers�o Melhorada
 
 P(6,1)=0.3; P(6,11)=0.45; P(6,15)=0.25; %N� 6
 P(15,20)=0.25; P(15,5)=0.40; P(15,6)=0.35; %N� 15
 P(1,6)=0.20; P(1,20)=0.25; P(1,7)=0.25; P(1,13)=0.30; %N� 1
 P(3,12)=0.4; P(3,19)=0.3; P(3,16)=0.3; %N� 3
 P(12,3)=0.30; P(12,10)=0.25; P(12,8)=0.45; %N� 12
 P(4,19)=0.4; P(4,13)=0.15; P(4,2)=0.45; %N� 4
 P(20,15)=0.35; P(20,1)=0.15; P(20,7)=0.15; P(20,14)=0.35; %N� 20
 P(13,1)=0.25; P(13,2)=0.35; P(13,4)=0.15; P(13,19)=0.25; %N� 13
 P(16,3)=0.25; P(16,7)=0.15; P(16,18)=0.6; %N� 16
 P(10,17)=0.45; P(10,12)=0.35; P(10,9)=0.2; %N� 10
 
 erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,d,P,0,'a','a',0);

 erro_med=erro./n_runs;
 t=1:n_passos;
 
 figure(2)
 plot(t, erro_med);
 hold on
 
 % Vers�o Piorada
 
 clear P
 
 load('MarkovChain.mat');

 P(1,6)=0.8; P(1,20)=0.1; P(1,7)=0.1;
 P(6,1)=0.1; P(6,11)=0.45; P(6,15)=0.45;
 
 erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,d,P,0,'a','a',0);

 erro_med=erro./n_runs;
 t=1:n_passos;
 
 figure(2)
 plot(t, erro_med);
 xlabel('t [s]'); ylabel('Erro');
 title('Erro de estimativa de posi��o da fonte');
 legend('Melhor circula��o','Pior circula��o');
%%
% A partir dos gr�ficos obtidos, pode-se verificar, que quando se inicia o
% token numa posi��o central, por exemplo �ncora 7, o erro diminui
% exponencialmente nos primeiros instantes. No entanto, neste caso, converge para
% um valor diferente de 0, ou seja a estima��o da localiza��o da fonte n�o
% � exata. 
% 
% Do mesmo modo, determinou-se a situa��o em que se inicia o token no estado 15, 
% que � uma �ncora pertencente a um subconjunto de estados 
% que tende a reter o mesmo e que se encontra longe da fonte. Como era de
% esperar neste caso, o erro leva mais tempo a aproximar-se da fonte, no
% entanto, converge aproximadamente para o mesmo valor que se verificou
% na situa��o anterior.
%
% Por outro lado, quando se inicia o token na �ncora 10, 
% que tamb�m pertence a 
% uma zona que tende a reter o mesmo, no entanto bastante mais
% pr�xima da fonte, este caso, relativamente �s duas outras situa��es, � o mais
% r�pido a convergir e converge para um erro=0.
%
% A partir da segunda figura, conclu�mos que as variantes do grafo da
% al�nea 2.d) influenciam a evolu��o do erro. Tal como era esperado, o grafo
% melhorado no sentido da equipa converge mais rapidamente para a posi��o 
% exata da fonte, pois o token circula de forma eficaz por todo o
% grafo. 
%
% Para o grafo onde se dificulta a circula��o do token, verifica-se que
% este converge para um erro bastante superior a 0, uma vez que n�o
% consegue chegar a todas as antenas, logo dificulta a estimativa exata da
% posi��o da fonte.

%% Quest�o 3.c)
% Nesta al�nea, a fonte movimenta-se em simult�neo com a transi��o do token.
% Desta forma, introduz-se um factor de esquecimente 0< $\lambda$ $\leq$ 1 no
% algoritmo RLS, para que na fun��o de custo se d� mais peso aos �ltimos
% termos do que os primeiros.
%%

clear
 close all
 
 load('MarkovChain.mat')
 n_passos=1200;
 n_runs=300;

 Po=100; %Pot�ncia da fonte
 variancia=10^-2;
 lambda_set=[1 0.3 0.85];
 
 a=[nodePos(:,2),nodePos(:,3)]';
 D=squareform(pdist([sourcePos' zeros(size(sourcePos')) a]'));
 d=D(1,3:end); %Fonte-�ncora dist�ncia

 for i=1:length(lambda_set)
     
     erro=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,...
                          sourcePos,d,P,0,'a','fonte_movimento',lambda_set(i)); 
     erro_med=erro./n_runs;
     t=1:n_passos;
     
     figure(1)
     plot(t, erro_med);
     hold on
 end


figure(1)
xlabel('t [s]'); ylabel('Erro'); 
legend('\lambda=1','\lambda=0.3','\lambda=0.85');
title('Evolu��o do erro para uma fonte em movimento')

%%
% * Para $\lambda$=1, pode-se observar que o erro diminui exponencialmente
% numa primeira parte, chegando a 0 por uns instantes. No entanto, a partir
% de um certo n�mero de transi��es come�a a aumentar, devido
% ao factor de esquecimento dar o mesmo peso a todas as medi��es. Com a
% fonte em movimento, � necess�rio introduzir um $\lambda$ de forma a
% dar-se mais peso �s medi��es mais recentes em compara��o com as antigas.
%
% * Para um valor baixo de $\lambda$, por exemplo 0.3, o algoritmo d�
% import�ncia apenas � �ltima medi��o, praticamente ignorando as medi��es
% anteriores. Como se pode verificar, esta estrat�gia tamb�m n�o conduz a uma
% boa estimativa da fonte. � necess�rio encontrar um equil�brio entre o peso
% das medi��es actuais e as anteriores.
%
% * De seguida, testaram-se v�rios valores de $\lambda$ de forma a
% encontrar um valor que melhor estimasse a posi��o da fonte para todas as
% transi��es, obtendo-se $\lambda$=0.85.
%
% Por fim, � necess�rio ter em aten��o que o erro pode ser maior ou menor
% tendo em conta a traject�ria da fonte e se esta se aproxima ou afasta
% das antenas.
