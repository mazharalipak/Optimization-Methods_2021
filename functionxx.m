function fun =functionxx (V)

X = V(1,:);
Y = V(2,:);

fun  = sin(Y).* exp(1- cos(X)).^2 + cos(X).* exp(1-sin(Y)).^2 + (X-Y).^2;

%X.^2 - 4*X + Y.^2 - Y - X.*Y;

end