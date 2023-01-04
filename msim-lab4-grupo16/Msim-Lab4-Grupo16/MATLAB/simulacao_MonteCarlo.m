function estados_MC= simulacao_MonteCarlo(P,n_runs,n_passos,n_estados,mode,estado_inicial)
%Entradas: P - matriz probabilidade de transição entre estados
% n_estados - número total de estados 
% mode - escolha entre uma inicialização sempre no mesmo estado_inicial para
% cada run (modo 'c') ou aletória (uniforme)
%Saídas: matriz que guarda o número total de vezes que o token passou por
%        cada antena, por passo.

estados_MC=zeros(n_passos, n_estados);

for i=1:n_runs
    if(mode=='c') %verificação do modo
        estado_atual=estado_inicial;
    else 
        estado_atual=randi(n_estados,1);
    end
    
    for j=1:n_passos
        %incrementar a matriz para o passo j e a âncora onde se encontra
        estados_MC(j, estado_atual)=estados_MC(j, estado_atual)+1; 
        %Escolha do estado seguinte
        estado_atual=find(cumsum(P(estado_atual,:)) > rand,1,'first');
    end 
    hh=waitbar(i/n_runs);
end

close(hh);
end