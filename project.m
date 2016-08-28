function g = project (x,eta)
T = eta*x.V;
g.M = x.U'*T; 
g.Up = T - x.U*(x.U'*T);
g.Vp = eta'*x.U; 
g.Vp = g.Vp - x.V*(x.V'*g.Vp);
end


