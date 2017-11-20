function dx = dRK4(x, t, h)

f = @L96;
df = @dL96;

dx = eye(length(x));

x = RK4(t, x, h);

for i = 1:length(t)-1
    k1 = f(x(:,i));
    dk1 = df(x(:,i))*dx;
    k2 = f(x(:,i) + 1/2*h*k1);
    dk2 = df(x(:,i) + 1/2*h*k1) * (dx + 1/2*h*dk1);
    k3 = f(x(:,i) + 1/2*h*k2);
    dk3 = df(x(:,i) + 1/2*h*k2) * (dx + 1/2*h*dk2);
    k4 = f(x(:,i) + h*k3);
    dk4 = df(x(:,i) + h*k3) * (dx + h*dk3);
    
    dx = dx + 1/3*h*(1/2*dk1 + dk2 + dk3 + 1/2*dk4);
end

end