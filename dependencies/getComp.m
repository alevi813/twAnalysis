function comp = getComp 

[~, host] = system('hostname');

if strcmp(host(1), 'd') %dhcp-129-116-178-237.cps.utexas.edu
    comp = 'desktop';
else
    comp = 'laptop';
end