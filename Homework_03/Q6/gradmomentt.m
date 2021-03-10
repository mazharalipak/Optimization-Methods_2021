%% Q8 ......... Steepest descent method ..........

function [Iter, Tole, V_k , alpha, XX, del_XX] = gradmomentt(V)

%% Initialization ......................

nx    = 5;
A     = hilb(nx);

%% parametrs ...................

beta  = 0.5;
   LX = 2*norm(A,'fro');
alpha = ( 1+beta )/LX;

%%

f = functionxx(A);                                                             % Original function .......................................

%%

V_k      = V;

%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-4)  
 
    V_kprev = V_k;    
    V_kk    = V_k - alpha*grada( V_k, A ) + beta*( V_k - V_kprev );
    V_k     = V_kk;    

%% storing variables .........................

XX(:,Iter) = V;                                                            

del_X          = V_k - V_kprev;
del_XX (Iter)  = norm( del_X ); 

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
title ('Q6: Heavy Ball Method','Interpreter','Latex','fontsize',14);
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

