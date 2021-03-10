%% Q8 ......... Steepest descent method ..........

function [Iter, Tole, V_k , alpha, XX, del_XX] = gradfast(V)

%% Initialization ......................


alpha = 1.0;                                                               % Setting stepsizes .................................................
nx    = 5;
A     = hilb(nx);

%%

f = functionxx(A);                                                             % Original function .......................................

%%

theta_k = 1;
y_k     = V;
V_k     = V;


%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-4)  

    LX = 2*norm(A,'fro');
 alpha = 1/LX;
 
 
 %%
 
 theta_prev  = theta_k;
 V_kprev     = V_k;
 V_k         = y_k - alpha*grada( y_k, A );
 theta_k     = ( 1 + sqrt( 1 + 4*theta_k^2 ) )/2;
 y_k         = V_k + ( ( theta_prev - 1 )/theta_k )*( V_k - V_kprev );

%% storing variables .........................

XX(:,Iter) = V;                                                            

del_X          = V_k - V_kprev;
del_XX (Iter)  = norm(del_X); 

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole = norm(del_X);   
    counter = counter + 1;
    if counter == 1000
        break;
    end 
end

%% plotting ...............


semilogy(del_XX,'b','LineWidth',1.0)
title ('Q6: Fast Gradient Method','Interpreter','Latex','fontsize',14);
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


