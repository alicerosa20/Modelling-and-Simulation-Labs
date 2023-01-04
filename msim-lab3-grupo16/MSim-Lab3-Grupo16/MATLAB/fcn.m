function d2teta = fcn(teta,dteta, T, m, l, M, L, g, k, beta)
%#codegen
J=(m*l^2)+(M/3)*L^2;
d2teta =(1/J)*(-k*teta*(1+(teta^2/100))-beta*dteta+g*sin(teta)*(m*l+(1/2)*M*L)+T);
