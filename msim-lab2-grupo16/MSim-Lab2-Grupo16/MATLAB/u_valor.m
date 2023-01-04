function u=u_valor(y,dy,yl)

k1=1/yl;
k2=sqrt(2*k1);

 if abs(-y)<=yl
     f=(k1/k2)*(-y);
 elseif abs(-y)>yl
     f=sign(-y)*(sqrt(2*abs(-y))-(1/k2));
 end
 
 u=sign(f-dy);
 
end