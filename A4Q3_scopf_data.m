% Modified IEEE 9 bus test system WECC representation
Sbase = 100; % MVA base

% Network information
ifrom = [1 4 5 3 6 7 8 8 9]';
ito = [4 5 6 6 7 8 2 9 4]';
x = [0.0576 0.092 0.17 0.0586 0.1008 0.072 0.0625 0.161 0.085]';
is = 1;
fmax = [250 250 150 300 150 250 250 250 105]'; % Max branch MW flows

% Nodal demands
d = [0 0 0 0 200 0 100 0 225]';

% Generator information
co = [800 400 400]';
a = [20 23 56]';
b = [.030 .035 .04]';
gmax = [250 300 270]';
gmin = [10 10 10]';
ngen = [1 2 3]';
dp = 1; % for LMP calculation