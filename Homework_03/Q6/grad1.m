%% Q8 ......... Steepest descent method ..........

function [Iter, Tole, V , alpha, XX, del_XX] = grad1(V)

%% Initialization ......................


alpha = 1.0;                                                               % Setting stepsizes .................................................
nx    = 5;
A     = hilb(nx);


f = functionxx(A);                                                             % Original function .......................................

%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;

while ( Tole > 1e-6 )  

X_ini = V;                                                                 % vector of variables at previous step .......................                                        

grad_fun = grada( V, A );                                                     % Gradient of the function ...........................................

%% Constant alpha

   LX = 2*norm(A,'fro');
   alpha = 1/LX;

%% Steepest descent 

%   alpha = argmin(grad_fun,A);

%% Armijo inexact line search ..............

%   alpha = armijo( V, grad_fun, f);                                              % Finding optimal alpha ............................

%%

V = V - alpha*( grad_fun );                                                % Updating the variables ........................

%% storing variables .........................

XX(:,Iter) = V;                                                            

del_X          = V - X_ini;
del_XX (Iter)  = norm(del_X); 

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole = norm( del_X );   
    counter = counter + 1;
    if counter == 1000
        break;
    end 
end

%% plotting ...............



semilogy(del_XX,'b','LineWidth',1.0)
title ('Q6: GD with exact line search (steepest descent)','Interpreter','Latex','fontsize',14);
xlabel('Iterations','Interpreter',' Latex','fontsize',14);
ylabel('$||x^{k+1} - x_{k}||$','Interpreter',' Latex','fontsize',14);    
xlim([-Inf Inf])
ylim([-Inf Inf])



%% Function ..............................

function f = functionxx(A)

f = @(X) X'*A*X;


%% Gradient of function f .......................

function [grad_fun] = grada( X, A )

grad_fun = 2*A*X;

%% Finding optimum alpha by Armijo line search ..................................

function alpha = armijo( X , g, f)

alpha   = 1;
fac     = 0.9;                                                        % < 1 reduction factor of alpha
c_1     = 0.1;
        
while f ( X - alpha*g ) > ( f( X ) - c_1*alpha*(g'*g) )
    
    alpha = fac*alpha;
    
    if alpha < 10*eps
        error('Error in Line search - alpha close to working precision');
    end  
    
end


function alpha = argmin(g,A)

alpha = (g'*g)/ ( 2*(g'*A*g) );


