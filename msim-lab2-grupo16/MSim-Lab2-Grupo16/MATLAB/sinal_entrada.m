function [u,t]=sinal_entrada(T,alpha,beta,U1,U2,n1,n2)

    T1=T/(1+alpha);
    T2=alpha*T1;
    t1=linspace(-(1/2)-(beta/2),(1/2)+(beta/2),n1);
    t2=linspace(-(1/2)-(beta/2),(1/2)+(beta/2),n2);
    
    u1=-U1*imp_prot(t1,beta); %Mudar a amplitude do sinal
    t_1=(t1*(T1/(1+beta)))+(T1/2); %Expansão e deslocamento do sinal do tempo
    
    u2=U2*imp_prot(t2,beta); 
    t_2=(t2*(T2/(1+beta))); %Expansão do sinal do tempo
    dif=t_1(n1)-t_2(1); %Distância entre o final do 1º impulso e inicio do segundo
    t__2=t_2+dif; %Deslocamento do segundo impulso no tempo
    
    u=[u1 u2]; %Concatenação dos dois impulsos 
    t=[t_1 t__2];
end