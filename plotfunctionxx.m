%% Plot function ................

function [XX, YY, funC]  = plotfunctionxx (Radius)

RR = (0:0.1:Radius)';
theta = (2*pi*(0:100))/100;

XX = -Radius + RR*cos(theta);
YY = -Radius + RR*sin(theta);

funC = sin(YY).* exp(1- cos(XX)).^2 + cos(XX).* exp(1-sin(YY)).^2 + (XX-YY).^2;

end
