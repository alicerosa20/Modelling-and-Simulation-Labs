function media_BPM=calc_media_BPM(t,locs,ind_inicial)
k=1;

for i=ind_inicial:(length(locs)-1)
    f_osc=1/(t(locs(i+1))-t(locs(i)));
    BPM(k)=60*(2*f_osc);
    k=k+1;
end

media_BPM=mean(BPM);

end