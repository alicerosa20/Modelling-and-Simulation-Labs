function x=sourcePos_estimativa(M,Po,variancia,Vp,a,d)

% Inicialização das matrizes A, b e Pi
A=zeros(M,4);
b=zeros(M,1);
Pi=zeros(M,1);

token_anc=round(Vp.*M); % Número de vezes que o token passa por cada âncora

% calcula-se a diferença entre o número de total de medições e o número de
% medições definido. Se a diferença for negativa, soma-se o módulo ao
% número de medições da última âncora. Se for positiva, substrai-se.

dif=sum(token_anc)-M;

if dif<0
    token_anc(20)=token_anc(20)+abs(dif);
elseif dif>0
    token_anc(20)=token_anc(20)-dif;
end

n=randn(M,1)*sqrt(variancia); % distribuição gaussiana do ruído
j=1;

for i=1:20
    contagem=token_anc(i); % Contagem do número total de medições para cada âncora 
    while contagem>0 % Contagem decrescente 
        Pi(j)=(Po/(d(i))^2)*exp(n(j)); %Potência na observação j
        %construção das matrizes A e b
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

x=[w(1) w(2)]; %estimativa posição da fonte

end