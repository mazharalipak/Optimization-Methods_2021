% Nedler Mead Method .........................................
 
function [Iter, Tole, V] = Nelder_new(G, B, W, alpha, gamma,beta, Radius)
%% Iteration counter ........................

% alpha = 1;
% gamma = 0.5;
% beta  = 2;

V = [B G W] ;

%%

Tole = 1;  
Iter = 1;
counter = 0;

tic;
while (Tole > 1e-3)
%% Tolerance criteria....
        
scale =  norm(B-G) + norm(W-B);   
        
%% Evaluating function ..............

fun = functionxx (V);

%% Sorting the simplex with Best (B) Good (G) and worst (W) ......................... 

[~ , sorb] = sort(fun);
V = (V(:,sorb));
B = V(:,1);
G = V(:,2);
W = V(:,3) ;

%% Reflection ............

M =  (B + G)./2;                             % mid point .............

R = M + alpha*(M - W);                              % Reflection .................
   
%% Logical operation.....................................

% Expansion ...............

if functionxx(R) < functionxx(W) 
    V(:,3) = R; 
E = R + gamma*(R-M);                               % Expansion point 
if functionxx(E) < functionxx(R)
    V(:,3) = E;
else
    V(:,3) = R;
end
    
% contraction ......................

elseif functionxx(W) < functionxx(R)
C1 = (W + M)/2;
C2 =  C1 + beta*(M-C1);
if functionxx(C1) < functionxx(C2) && functionxx(C1)< functionxx(W)
    V(:,3) = C1;
elseif  functionxx(C2) < functionxx(C1) && functionxx(C2)< functionxx(W)
    V(:,3) = C2;
else
    S = (B+W)/2;
    V(:,3) = S;
    V(:,2) = M;
end
end
 
%% Storing data .................................

VH(1,Iter) = {cat(2,V,V(:,1))};

%% Iteration counter ...............................

Iter = Iter + 1;
    Tole=(abs(scale));   
    counter = counter + 1;
    if counter == 100
        break;
    end 
end

LL = cell2mat(VH);

%% Plotting domain and Iterations 

[XX, YY, funC]  = plotfunctionxx (Radius);

hold on

contour(XX,YY,funC)

plot( (LL(1,:)), (LL(2,:)),'Color','r','LineWidth',0.5,'MarkerSize',5,'MarkerEdgeColor','b')

hold on

xlim([-Inf Inf])
ylim([-Inf Inf])
  
title ('Question 3: Nedler Mead Iterations','Interpreter','Latex','fontsize',14);
xlabel('$x$','Interpreter',' Latex','fontsize',14);
ylabel('$y$','Interpreter',' Latex','fontsize',14);    
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.4 4.4])
print -djpeg filename.jpg -r600
end
