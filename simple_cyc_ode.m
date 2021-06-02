function dcyc = simple_cyc_ode(v,adj_beta,a,b,c,d,e,delta)

x = v(1:2:2*size(adj_beta,1)); %set the values for x in each nurse cell
y = v(2:2:2*size(adj_beta,1)); %set the values for y in each nurse cell


%solve for z at each step under the quasi-steady state approx.
z = inv(eye(size(adj_beta,1)) - delta*adj_beta)*[sum(x); zeros(size(adj_beta,1)-1,1)];

%update x and y for each cell
dx = (a+x.^2)./(1+x.^2)./(1+y)./(1+z)-b.*x;
dy = c - d.*y./(1+e.*(x.^2));

%no change in the oocyte (no limit cycle in this cell)
dx(1) = 0;
dy(1) = 0;

%update x and y for the next iteration
dcyc(1:2:2*size(adj_beta,1),1) = dx;
dcyc(2:2:2*size(adj_beta,1),1) = dy;
end