function small_example()
clc;
clear;
randn('state',0); rand('state',0);

m = 1000; n = 1000; k = 40; OS = 3;

LRGeom.m=m;
LRGeom.n=n;
LRGeom.r=k;
LRGeom.OS=OS;

% random factors
L = randn(m, k); 
R = randn(n, k); 
dof = k*(m+n-k);

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
prob = make_prob(L,R,Omega,k); % <- you can choose another rank here

options = default_opts();
x0 = make_start_x(prob);

t=tic;
[Xcg,hist_Fast] = LRGeomFOA(prob,options,x0);
out_time = toc(t)

