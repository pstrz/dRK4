function x = RK4(t, x, h)
% Runge-Kutte for L96

f = @L96;

for i = 1:length(t)-1
    k1 = f(x(:,i));
    k2 = f(x(:,i) + 1/2*h*k1);
    k3 = f(x(:,i) + 1/2*h*k2);
    k4 = f(x(:,i) + h*k3);
    
    x(:,i+1) = x(:,i) + 1/3*h*(1/2*k1 + k2 + k3 + 1/2*k4);
end

end