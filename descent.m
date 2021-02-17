%% Coordinate descent method ..........

function [Iter, Tole, V] = descent(V)

%% Initialization ......................

% Setting stepsizes
gamma = 0.008;
alpha = 0.05;

num = size(V,1); 

h = zeros(num,1);
%% Defining Iteration structure ........................................

Tole = 1;  
Iter = 1;
counter = 0;
tic;
while (Tole > 1e-5)  

F_not = functionxx(V);                                                     % Evaluating function at previous step .........................

X_ini = V;                                                                 % vector of variables at previous step .......................

% choosing the descent direction, x_{1}, x_{2},......

jj = mod(Iter,num);    
if (jj == 1)        
    h = [1;0];    
elseif (jj == 0)        
    h = [0;1];    
end

%% iteration process .....................................

fun = functionxx ( V + alpha*h );

gg = ( 1/alpha )*(fun'- F_not').*h;                                        

V = V - gamma*( gg);                                                       % Updating the variables ........................

F_final = functionxx(V);

scale = X_ini -V;

%% storing variables .........................

XX(:,Iter) = V;                                                            

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole=max(abs(scale));   
    counter = counter + 1;
    if counter == 1000
        break;
    end 
end


%% plotting ...............

[XX1, YY1, funC]  = plotfunctionxx (5);

hold on

contour(XX1,YY1,funC)

plot(XX(1,:),XX(2,:),'r.--')

title ('Question 4: Coordinate descent Iterations','Interpreter','Latex','fontsize',14);
xlabel('$x$','Interpreter',' Latex','fontsize',14);
ylabel('$y$','Interpreter',' Latex','fontsize',14);    
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.4 4.4])
print -djpeg filename.jpg -r600
 
end