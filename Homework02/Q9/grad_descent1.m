%% Q9 for n = 10  ......... Steepest descent method ..........

function [Iter, Tolex, V , alpha, XX] = grad_descent1(V)
%% Initialization ......................

alpha = 10*1e-3;                                                           % Setting stepsizes .................................................

f = functionxx;                                                            % Original function .......................................

%% Defining Iteration structure ........................................

Tolex = 1;  
Iter = 1;
counter = 0;

while (Tolex > 1e-6)  

X_ini = V;                                                                 % vector of variables at previous step .......................                                        

grad_fun = grada( V );                                                     % Gradient of the function ...........................................

alpha = opt_ff (V,grad_fun,f);                                             % Finding optimal alpha ............................

V = V - alpha*( grad_fun );                                                % Updating the variables ........................

%% storing variables .........................

XX(:,Iter) = V;                                                            

%%

g_tole = ( norm (V - X_ini ) );
g_tolee (Iter) = g_tole;
%% Iteration counter ...............................

Iter = Iter + 1;
    Tolex = g_tole ;   
    counter = counter + 1;
    if counter == 100
        break;
    end 
end

%% plotting ...............

%% plotting ...............


semilogy(g_tolee,'LineWidth',1.5,'MarkerSize',5,'MarkerEdgeColor','b');
xlim([-Inf Inf])
ylim([-Inf Inf])


%% Function ..............................

function f =functionxx

%% n = 10 ........................ 

 f = @(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) 0.25*( x1 - 1).^2 + ( 2*x1.^2 - x2 -1).^2 + ( 2*x2.^2 - x3 -1).^2 +...
( 2*x3.^2 - x4 -1).^2 + ( 2*x4.^2 - x5 -1).^2 +...
    ( 2*x5.^2 - x6 -1).^2 + ( 2*x6.^2 - x7 -1).^2 +...
    ( 2*x7.^2 - x8 -1).^2 + ( 2*x8.^2 - x9 -1).^2 + ( 2*x9.^2 - x10 -1).^2;

%% Gradient of function f .......................

function [grad_fun] = grada(X)

x1 = X(1,1); x2 = X(2,1); x3 = X(3,1); x4 = X(4,1); x5 = X(5,1);
x6 = X(6,1); x7 = X(7,1); x8 = X(8,1); x9 = X(9,1); x10 = X(10,1);

% gradient ..................

grad_fun = [ x1/2 - 8*x1*(- 2*x1^2 + x2 + 1) - 1/2;
    - 4*x9^2 + 2*x10 + 2;
    - 4*x1^2 + 2*x2 - 8*x2*(- 2*x2^2 + x3 + 1) + 2;
    - 4*x2^2 + 2*x3 - 8*x3*(- 2*x3^2 + x4 + 1) + 2;
    - 4*x3^2 + 2*x4 - 8*x4*(- 2*x4^2 + x5 + 1) + 2;
    - 4*x4^2 + 2*x5 - 8*x5*(- 2*x5^2 + x6 + 1) + 2;
    - 4*x5^2 + 2*x6 - 8*x6*(- 2*x6^2 + x7 + 1) + 2;
    - 4*x6^2 + 2*x7 - 8*x7*(- 2*x7^2 + x8 + 1) + 2;
    - 4*x7^2 + 2*x8 - 8*x8*(- 2*x8^2 + x9 + 1) + 2;
    - 4*x8^2 + 2*x9 - 8*x9*(- 2*x9^2 + x10 + 1) + 2];

%% Finding optimum alpha by solving min f( x - \alpha \nabla f )

function [alpha] = opt_ff (V,grad_fun,f)
syms aa

xx_new = V - aa*grad_fun;

Lx  = simplify( jacobian( (f( xx_new(1),xx_new(2),xx_new(3),xx_new(4),xx_new(5),xx_new(6),...
    xx_new(7),xx_new(8),xx_new(9),xx_new(10) )),aa ));
Cx = fliplr ( coeffs(Lx) );

alpha = min(real(double(roots(Cx))));


