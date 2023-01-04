function [erro,x_pos,x_est]=erro_MonteCarlo(n_passos,n_runs,Po,variancia,nodePos,sourcePos,...
                              d,P,estado_inicial,mode,mode_mov,lambda)

erro=zeros(n_passos,1);

for i=1:n_runs
    A=zeros(n_passos,4);
    b=zeros(n_passos,1);
    n=randn(n_passos,1)*sqrt(variancia);
    x = sourcePos';
    % Escolha do modo de uma condição inicial defenida ou random uniforme
    if(mode=='c')
        ancora=estado_inicial;
    else 
        ancora=randi(20,1);
    end
    % Escolha do modo da fonte em movimento,'fonte_movimento', ou fonte
    % parada e dependendo disso aplicar ou não o valor de Lambda
    if strcmp(mode_mov,'fonte_movimento')
           rlsPar=struct('lam',lambda);
    else 
           rlsPar=struct('lam',1);
    end
    
    for j= 1:n_passos
        %Verificar se a fonte está parada ou em movimento e dependendo
        %disso mudar as coordenadas da fonte em cada passo
        if strcmp(mode_mov,'fonte_movimento')
           x=x+[-0.02;0.02];
        end
        Pi=(Po/(d(ancora))^2)*exp(n(j));
        ai = nodePos(ancora,2:3)';
        
        %construção matriz A e b
        A(j,1)=-2*Pi*ai(1);
        A(j,2)=-2*Pi*ai(2);
        A(j,3)=-1;
        A(j,4)=Pi;
        b(j)=-Pi*(norm(ai))^2;
        
        [~,w,rlsPar]=qrrls(A(j,:),b(j),rlsPar);
        erro(j)=erro(j)+norm(x-w(1:2));
        
       ancora = find(cumsum(P(ancora, :))>rand,1,'first');
    end
    hh=waitbar(i/n_runs);
end
    close(hh);
end