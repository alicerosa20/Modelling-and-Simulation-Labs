function x=sourcePos_estimativa(M,Po,variancia,Vp,a,d)

% Inicializa��o das matrizes A, b e Pi
A=zeros(M,4);
b=zeros(M,1);
Pi=zeros(M,1);

token_anc=round(Vp.*M); % N�mero de vezes que o token passa por cada �ncora

% calcula-se a diferen�a entre o n�mero de total de medi��es e o n�mero de
% medi��es definido. Se a diferen�a for negativa, soma-se o m�dulo ao
% n�mero de medi��es da �ltima �ncora. Se for positiva, substrai-se.

dif=sum(token_anc)-M;

if dif<0
    token_anc(20)=token_anc(20)+abs(dif);
elseif dif>0
    token_anc(20)=token_anc(20)-dif;
end

n=randn(M,1)*sqrt(variancia); % distribui��o gaussiana do ru�do
j=1;

for i=1:20
    contagem=token_anc(i); % Contagem do n�mero total de medi��es para cada �ncora 
    while contagem>0 % Contagem decrescente 
        Pi(j)=(Po/(d(i))^2)*exp(n(j)); %Pot�ncia na observa��o j
        %constru��o das matrizes A e b
        A(j,1)=-2*Pi(j)*a(1,i);
        A(j,2)=-2*Pi(j)*a(2,i);
        A(j,3)=-1;
        A(j,4)=Pi(j);
        b(j)=-Pi(j)*(norm(a(1:2,i)))^2;
        
        contagem=contagem-1;
        j=j+1;
    end

end

rls_Par=struct('lam',1);
[~,w,~]=qrrls(A,b,rls_Par); % Retorna o vector de coeficientes do filtro transversal convencional, w

x=[w(1) w(2)]; %estimativa posi��o da fonte

end