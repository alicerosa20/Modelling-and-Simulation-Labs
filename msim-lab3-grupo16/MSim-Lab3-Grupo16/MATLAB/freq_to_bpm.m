function BPM=freq_to_bpm(wn)

f_osc=wn/(2*pi);
f_bat=2*f_osc;
BPM=60*f_bat;

end