V = [-25:25:25];
X = [-0.5:0.1:1.5];
for i = 1:length(V)
  v = V(i);
  for j = 1:length(X)
    x = X(j);
    dx(i,j) = mk_cr (x,v);
  end
end
figure(1);
plot(X,dx);
grid;
legend(num2str(V'))
xlabel("x");
ylabel("dx/dt");

## Test steady state.
[dx,x_0] = mk_cr (0,0)
dx = mk_cr (x_0,0)

## dx against v at steady state corresp v=0.
v = [-50:50];
for i = 1:length(v)
  dx_v(i) = mk_cr (x_0,v(i));
end

figure(2);
plot(v,dx_v);
grid;
xlabel("v");
ylabel("dx/dt");
