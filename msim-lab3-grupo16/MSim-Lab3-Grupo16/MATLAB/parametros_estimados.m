function [BPM_obt,ind_r,l,m]=parametros_estimados(BPM,valor_esperado,m_set,l_set)

[r,c] = size(BPM);
erro_min=50;

for i=1:r
    for j=1:c
        if isnan(BPM(i,j)) %Se o valor de BPM não é válido continua
            continue;
        end
         %Erro entre a frequência pretendida e a medida para o par (m,l)
        erro=abs(BPM(i,j)-valor_esperado);
        
        %Guardar os indíces que levam ao erro mínimo
        if erro<erro_min
            ind_r=i;
            ind_c=j;
            BPM_obt=BPM(i,j);
            erro_min=erro;
        end
    end
end

m=m_set(ind_r);
l=l_set(ind_c);

end