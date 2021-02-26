%% Q8 ......... Steepest descent method ..........

function [Iter, Tole, V , alpha, XX] = grad_deent(V)

%% Initialization ......................


alpha = 10*1e-3;                                                           % Setting stepsizes .................................................

f = functionxx;                                                             % Original function .......................................

%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-6)  

X_ini = V;                                                                 % vector of variables at previous step .......................                                        

grad_fun = grada( V );                                                     % Gradient of the function ...........................................

alpha = opt_ff (V,grad_fun,f);                                             % Finding optimal alpha ............................

V = V - alpha*( grad_fun );                                                % Updating the variables ........................


%% storing variables .........................

XX(:,Iter) = V;                                                            

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole = norm(grad_fun);   
    counter = counter + 1;
    if counter == 100
        break;
    end 
end

%% plotting ...............

figure(1); clf; ezcontour(f,[-10 10 -10 10]);

title ('$x^2 + xy + 10y^2 -22y -5x$','Interpreter','Latex','fontsize',14);
xlabel('$x$','Interpreter',' Latex','fontsize',14);
ylabel('$y$','Interpreter',' Latex','fontsize',14);    
axis equal;
hold on
plot( XX(1,:),XX(2,:),'ko-')


%% Function ..............................

function f =functionxx

f = @(x1,x2) x1.^2 + x1.*x2 + 10*x2.^2 - 22*x2 - 5*x1;


%% Gradient of function f .......................

function [grad_fun] = grada(X)

x1 = X(1,1);
x2 = X(2,1);

grad_fun = [ (2.*x1 + x2 - 5); 
    (x1 + 20.*x2 -22) ];

%% Finding optimum alpha by solving min f( x - \alpha \nabla f )

function [alpha] = opt_ff (V,grad_fun,f)
syms aa

xx_new = V - aa*grad_fun;

alpha = double( solve( jacobian( f(xx_new(1),xx_new(2) ) ) ) );
