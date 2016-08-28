function [x,histout,fail] = LRGeomFOA(prob, opts, x0)
 
%   Uses Armijo rule, polynomial linesearch on embedded submanifold of
%   fixed rank matrices.
%           
% Input: prob     = problem instance, see MAKE_PROB.
%        opts     = options, see DEFAULT_OPTS.           
%        x0       = starting guess.
%
% Output: x       = solution.
%         histout = iteration history. Each row of histout is
%                   [rel_norm(grad), rel_err_on_omega, relative_change, ...
%                        number of step length reductions, restarts]
%         fail    = flag for failure
%
% See also SMALL_EXAMPLE, MAKE_PROB, DEFAULT_OPTS. 


t_begin = tic();

%beta_type = 'F-R';
%beta_type = 'P-R';

fail = true;

norm_M_Omega = norm(prob.data);

% ORTH_VALUE = 0.1; % the search directions should be almost orthogonal

itc = 1; xc = x0;
fc = F(prob,xc);
gc = grad(prob,xc);
%ip_gc = ip(xc,gc,gc);
% fold = 2*fc; %reschg = -1; 
% first search-dir is steepest gradient
dir = scaleTxM(gc,-1);
% rel_grad = sqrt(ip_gc)/max(1,norm(xc.sigma));

ithist=zeros(opts.maxit,2);

k = length(xc.sigma);
yc = xc;
t = 1;
%xmat = xc.U*diag(xc.sigma)*xc.V';

for itc=1:opts.maxit

  tinit =2.5; %m =5000,n=5000,r =40,os=5,tinit=7.2;
 
    [xc_new,fcx_new,succ,numf,iarm,tinit] = armijo_search(prob, yc,fc,gc,dir, tinit);
 
   if itc ==1;
       yc_new =xc_new;
       fcy_new = fcx_new;
        newt = (1+sqrt(1+4*t^2))/2;
        gc_new = grad(prob,yc_new); 
        ip_gc_new = ip(yc_new,gc_new,gc_new);
         dir = scaleTxM(gc_new, -1);
   %       newxmat = xc_new.U*diag(xc_new.sigma)*xc_new.V';
   elseif itc ==2;
       yc_new =xc_new;
       fcy_new = fcx_new;
        newt = (1+sqrt(1+4*t^2))/2;
        gc_new = grad(prob,yc_new); 
        ip_gc_new = ip(yc_new,gc_new,gc_new);
         dir = scaleTxM(gc_new, -1);
          newxmat = xc_new.U*diag(xc_new.sigma)*xc_new.V';
   elseif itc>2
      
        newt = (1+sqrt(1+4*t^2))/2;
        beta = (t-1)/newt;
        newxmat = xc_new.U*diag(xc_new.sigma)*xc_new.V';
        
        newx_x = project(xc_new,newxmat-xmat);
        yc_new =  moveEIG(prob,xc_new,newx_x,beta);
        fcy_new = F(prob,yc_new);
        
        if fcy_new- fcx_new >0
          yc_new =xc_new;
          gc_new = grad(prob,yc_new);
          ip_gc_new = ip(xc_new,gc_new,gc_new);
          dir = scaleTxM(gc_new, -1);
          fcy_new =fcx_new;
        else
            gc_new = grad(prob,yc_new); 
             ip_gc_new = ip(yc_new,gc_new,gc_new);
               dir = scaleTxM(gc_new, -1);
        end
   end    
  % Test for convergence
  if fcy_new<1e-12
      fail = false;
      disp('f tol reached')
      break;
  end
  if sqrt(2*fc) < opts.abs_f_tol
    if opts.verbosity > 0
      disp('Abs f tol reached.')
    end
    fail = false;
    break;
  end
  if sqrt(2*fc)/norm_M_Omega < opts.rel_f_tol
    if opts.verbosity > 0; disp('Relative f tol reached.'); end
    fail = false;
    break;
  end
     
  if opts.stagnation_detection && itc > 10
    R1 = qr([xc_new.U*diag(xc_new.sigma) xc.U*diag(xc.sigma)],0);
    R1 = triu(R1(1:k,:));
    R2 = qr([xc_new.V -xc.V],0);
    R2 = triu(R2(1:k,:));
    
    %rel_change_x = norm(xc_new.U*diag(xc_new.sigma)*xc_new.V' - xc.U*diag(xc.sigma)*xc.V', 'fro') / ...
    %  norm(xc_new.U*diag(xc_new.sigma)*xc_new.V', 'fro');
    rel_change_x = norm(R1*R2','fro')/norm(xc_new.sigma,2);

    if rel_change_x < opts.rel_tol_change_x
      if opts.verbosity > 0; disp('Iteration stagnated for rel_tol_change_x.'); end
      fail = true;
      break;
    end
  end 
        
      

  % update _new to current
  t = newt;
  gc = gc_new;
  ip_gc = ip_gc_new;
  yc = yc_new;
  xc = xc_new;
  fc = fcy_new;
  rel_grad = sqrt(ip_gc)/max(1,norm(yc_new.sigma));
  if itc>1
  xmat = newxmat;
  end

 iter_end_time= toc(t_begin);
  ithist(itc,1) = iter_end_time; %rel_grad;
  ithist(itc,2) = fc;             %sqrt(2*fc)/norm_M_Omega;

end

x = xc_new;
ithist(itc,1) = iter_end_time;%rel_grad;
ithist(itc,2) = fc;%sqrt(2*fc)/norm_M_Omega;
histout=ithist(1:itc,:); 

