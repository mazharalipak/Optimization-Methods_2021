%% Q9 for n = 10  ......... Steepest descent method ..........

function [Iter, Tolex, V , alpha, XX] = grad_descent2(V)
%% Initialization ......................

alpha = 10*1e-3;                                                           % Setting stepsizes .................................................

f = functionxx;                                                            % Original function .......................................

%% Defining Iteration structure ........................................

Tolex = 1;  
Iter = 1;
counter = 0;

while (Tolex > 1e-5)  

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

f = @(x1,x2,x3) 0.25*( x1 - 1).^2 + ( 2*x1.^2 - x2 -1).^2 + ( 2*x2.^2 - x3 -1).^2 ;

%% Gradient of function f .......................

function [grad_fun] = grada(X)

x1 = X(1,1);
x2 = X(2,1);
x3 = X(3,1);


grad_fun = [ 0.5*( x1 -1 ) + 8*x1*( 2*x1^2 - x2 -1);
    -2*( 2*x1^2 -x2 - 1) + 8*x2*( 2*x2^2 - x3 -1);
    -2*( 2*x2^2 -x3 -1)];

%% Finding optimum alpha by solving min f( x - \alpha \nabla f )

function [alpha] = opt_ff (V,grad_fun,f)
syms aa

xx_new = V - aa*grad_fun;

alphax = double( solve( jacobian( (f(xx_new(1),xx_new(2),xx_new(3) )),aa ) ) );
alpha = min(abs(alphax));

