function out = deriv_eval(P,s,ord)

Pxs = fnder(P.Px,ord);
Pys = fnder(P.Py,ord);
Pzs = fnder(P.Pz,ord);
out = [ppval(Pxs,s);ppval(Pys,s);ppval(Pzs,s)];

end