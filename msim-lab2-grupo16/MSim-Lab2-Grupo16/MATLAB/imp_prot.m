function p_beta=imp_prot(t,beta)

ind= t<=-(beta/2)-(1/2);
part0=zeros(1,length(t(ind)));

ind= t>-(beta/2)-(1/2) & t<=-(1/2);
part1 = (1/2)*(((4/beta^2).*(t(ind)+(1/2))+(2/beta)).*(t(ind)+(1/2)+(beta/2)));

ind= t>-(1/2) & t<=(beta/2)-(1/2);
part2= 1-((1/2)*((-(4/beta^2).*(t(ind)+(1/2))+(2/beta)).*((beta/2)-(t(ind)+(1/2)))));

ind= t>(beta/2)-(1/2) & t<=-(beta/2)+(1/2);
part3=ones(1,length(t(ind)));

ind= t>-(beta/2)+(1/2) & t<=(1/2);
part4= 1-((1/2)*(((4/beta^2).*(t(ind)-(1/2))+(2/beta)).*(t(ind)-(1/2)+(beta/2))));

ind= t>(1/2) & t<=(beta/2)+(1/2);
part5=(1/2)*((-(4/beta^2).*(t(ind)-(1/2))+(2/beta)).*((beta/2)-(t(ind)-(1/2))));

ind= t>(beta/2)+(1/2);
part6=zeros(1,length(t(ind)));

p_beta=[part0 part1 part2 part3 part4 part5 part6];

end