clear
close all;
mex -setup C++;

mex MIMOalgo.cpp '-g -O3';
a = [1.1, 2.1, 3; 4, 5, 6; 7, 8, 9];
mex(a)
