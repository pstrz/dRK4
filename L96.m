function f = L96(x)

F = 8;
dim = length(x);
G = F * ones(dim,1);

f = (circshift(x,[-1,0]) - circshift(x,[2,0])).*circshift(x,[1,0]) - x + G;

end