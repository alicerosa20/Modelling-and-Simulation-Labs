function y=calc_erro_l2(x)

options = simset('SrcWorkspace','current');

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

y=abs(141-media_BPM);

end