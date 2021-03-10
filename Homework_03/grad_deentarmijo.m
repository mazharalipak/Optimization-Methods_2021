%% Q8 ......... Steepest descent method ..........

function [Iter, Tole, V , alpha, XX] = grad_deentarmijo(V)

% V - vector of system variables x_1, x_2 ................
% f - objective function .................
% grada_f  - gradient of objective function ........................
% alpha    - step size ...........................


%% Initialization ......................


alpha = 1.0;                                                               % Setting stepsizes .................................................

f = functionxx;                                                            % Original function .......................................

%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-3)  

X_ini = V;                                                                 % vector of variables at previous step .......................                                        

grad_f = grada( V );                                                     % Gradient of the function ...........................................

alpha = armijo( V, grad_f, f);                                              % Finding optimal alpha ............................

V = V - alpha*( grad_f );                                                % Updating the variables ........................


%% storing variables .........................

XX(:,Iter) = V;                                                            

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole = norm(grad_f);   
    counter = counter + 1;
    if counter == 100
        break;
    end 
end

%% plotting ...............

figure(1); clf; fcontour(f, [-0.5 1 -0.5 0.8]);

title ('$2x1^4 + 3x2^4 + 2x1^2 + 4x2^2 + x1x2 - 3x1 - 2x2$','Interpreter','Latex','fontsize',14);
xlabel('$x_1$','Interpreter',' Latex','fontsize',14);
ylabel('$x_2$','Interpreter',' Latex','fontsize',14);    
axis equal;
hold on
plot( XX(1,:),XX(2,:),'ko-')


%% Function ..............................

function f =functionxx

f = @(x1,x2) 2*x1^4 + 3*x2^4 + 2*x1^2 + 4*x2^2 + x1*x2 - 3*x1 - 2*x2;


%% Gradient of function f .......................

function [grad_fun] = grada(X)

x1 = X(1,1);
x2 = X(2,1);

grad_fun = [8*x1^3 + 4*x1 + x2 - 3
     12*x2^3 + 8*x2 + x1 - 2];

%% Finding optimum alpha by Armijo line search ..................................

function alpha = armijo( x, g, f)

alpha   = 1;
fac     = 0.9;                                                        % < 1 reduction factor of alpha
c_1     = 0.1;
        
while f ( x(1,1) - alpha*g(1,1), x(2,1) - alpha*g(2,1) ) > ( f ( x(1), x(2) ) - ( c_1*alpha*(g'*g) ) )
    
    alpha = fac*alpha;
    
    if alpha < 10*eps
        error('Error in Line search - alpha close to working precision');
    end  
    
end


