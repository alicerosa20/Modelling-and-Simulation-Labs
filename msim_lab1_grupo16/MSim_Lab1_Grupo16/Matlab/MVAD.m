function y=MVAD(x)

options = simset('SrcWorkspace','current');
load('presas.mat');

N1_inicial=4;
N2_inicial=x(1);
delta1=3.1;
delta2=-1.5;
alpha1=1.4;
alpha2=x(2);
tt=20;
x_max=-1;

sim('predador_presa',[],options);

 %Cálculo no máximo valor absoluto das diferenças
 for k=1:length(tr)
    
    x=abs(N1(k)-yr(k));
    
    if x>x_max
       x_max=x;
    end
    
 end

y=x_max;
 
end