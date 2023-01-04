function y=calc_erro_l1(x)
% Defenir novo workspace de simulação
options = simset('SrcWorkspace','current'); 
%Parâmetros
L=0.25;
M=0.1;
k=0.35;
g=9.8;
beta=0.001;
dteta_inicial=0;
teta_inicial=pi/4;
T=0;
m=0.0799;
l=x;

sim('sim10',40,options);

[~,locs] = findpeaks(teta);
ind_inicial=10;
media_BPM=calc_media_BPM(t,locs,ind_inicial);

y=abs(53-media_BPM);

end