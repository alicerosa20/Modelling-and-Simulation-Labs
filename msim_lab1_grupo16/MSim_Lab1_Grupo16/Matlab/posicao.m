 function teta=posicao(x,y,l)
 %Verificar se est� na origem 
 if x==0 && y==0 
     teta(1)=0;
     teta(2)=pi;
 else
     d=sqrt(x^2+y^2);
     sinal_x=x/abs(x); %guarda o sinal do x
     sinal_y=y/abs(y); %guarda o sinal do y
     
     if x==0 %Se o x=0 s� nos interessa o sinal do y
         teta(1)=(acos(d/(2*l))+acos(-y/d));
         teta(2)=acos((-y/l)-cos(teta(1)))*sinal_y;
     else
        %Se o sinal_x<0 queremos que o �ngulo seja negativo
        teta(1)=(acos(d/(2*l))+acos(-y/d))*sinal_x; 
        %Se o x<coordenada x da ponta da barra 1 ent�o teta_2 � negativo,
        %caso contr�rio teta_2 � positivo
        if x<l*sin(teta(1)) 
            sinal_x=-1;
        else
            sinal_x=1;
        end
        
        teta(2)=acos((-y/l)-cos(teta(1)))*sinal_x;
      end 
    end
 end
 