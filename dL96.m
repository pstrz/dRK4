function df = dL96(x)

X = diag(x);

df =  circshift(X, [-1,0]) - circshift(X, [2,3]) + circshift(X, [2,1]) - circshift(X, [-1,1]) - eye(length(x));

end