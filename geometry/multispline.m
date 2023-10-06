function [P,x] = multispline(x,s1,s2)

P.Px = spline(s1,x(1,:));
P.Py = spline(s1,x(2,:));
P.Pz = spline(s1,x(3,:));
x = [ppval(P.Px,s2);ppval(P.Py,s2);ppval(P.Pz,s2)];

end