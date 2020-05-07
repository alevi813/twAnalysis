% genvec=fitguass(lambda,t)
% 
% a = lambda(1);
% b = lambda(2);
% phi1 = lambda(3);
% 
% phi2 = abs(phi1 + pi);
% 
% sigma = lambda(4);
% dc = lambda(5);


function genvec=fitgauss(lambda,t) 

a = lambda(1);
b = lambda(2);
phi = lambda(3);
sigma = lambda(4);
dc = lambda(5); 

angle1 = angle(exp(1i.*(t - phi))); 
angle2 = angle(exp(1i.*(t - phi + pi)));

%this is wrong!?
genvec = dc + (a/(2*sigma^2)).*exp(-(angle1.^2)) + (b/(2*sigma^2)).*exp(-(angle2.^2));

% jake hack
%min + range*exp(-variance * (1 - cos(t - phi))

% phi1 = lambda(3);
% phi2 = abs(phi1 + pi);
% % 
% y1 = a * exp(-sigma*(1-sin(t - phi1)));
% y2 = b * exp(-sigma*(1-cos(t - phi2)));
% genvec = y1 + y2 + dc;


