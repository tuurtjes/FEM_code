%%% 
% cubic spline interpolation to find the shape of a EB element
%%%

x = [0 1];
y = [0 0.1];

dxdy  = [0 0.1];

spl = csape(x,[ 0 y 0.1],'complete');



fnplt(spl)