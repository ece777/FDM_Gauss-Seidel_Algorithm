clear all;


function curr = curr(n)
clc;
tol = 1e-6;
h = 180e-9;  # thickness
rho_bulk = 2.4318e-8; #including size-dependent resistivity
sigma = 1/rho_bulk;
%n = 7;
a = 15*sqrt(2)*1e-6;
b = 93*sqrt(2)*1e-6;
c = 78*sqrt(2)*1e-6;
rel = 1.9;

delta_y   = b/n;
nb        = round(b/delta_y); % Number of segments on ’b’
y1        = linspace(0,b,n+1); % Grid points along y-axis
na        = round(a/delta_y);
y2        = linspace(0,a,na+1); % still not sure

mc  = round(c/delta_y); % Number of segments on ’b’
z   = linspace(0,c,mc+1); % Grid points along y-axis


% Initialize potential and mask array
for k = 1:mc+1
  for j = 1:na+1+ (k-1)
    f(k,j) = 1;
    mask(k,j) = rel;
  end
end

% [mask(i,j) = 0 implies
% unchanged f(i,j)]
for k = 2:mc+1
  for j = 1:na+1+ (k-1)    
    f(k,j) = 0;
    end
end
k = 1;
  for j = 1:na+1+ (k-1)    
    mask(k,j) = 0;
  end
k = mc+1;
  for j = 1:na+1+ (k-1)    
    mask(k,j) = 0;
  end  

I_1_old = 0;
for iter = 1:5000  
  iter_index = iter   
  f = seidel(f, mask, na, mc);
  I_start = curr_dens_1(mc, na, f, delta_y, sigma, h);
  I_mid = curr_dens_mid(mc, na, f, delta_y, sigma, h); 
  I_end = curr_dens_2(mc, na, f, delta_y, sigma, h);
  I_sumavg = curr_dens_3(mc, na, f, delta_y, sigma, h)   
  R_start = 1 / I_start;
  R_mid = 1 / I_mid;
  R_end = 1 / I_end;
  R_sumavg = 1 / I_sumavg   # USE this!
  R   = (R_start + R_end)/2 ; #DO NOT use, not really acurate!
  
  E = E_field_1(mc, na, f, delta_y, sigma, n, iter) ;
  
  if ((abs((I_start - I_1_old)/I_start)) < tol)
    break    
  else
    I_1_old = I_start;   
  endif
  
endfor


#******************************************************************************
# To export to a file
result_matrix = ([f]);
dlmwrite (sprintf('test_grid_f_%dn_%diter.csv', n, iter), result_matrix); 
#
result_matrix = ([E]);
dlmwrite (sprintf('test_grid_E_%dn_%diter.csv', n, iter), result_matrix);
#
result_matrix = ([iter; I_start; I_end; R_start; R_end; R ]);
dlmwrite (sprintf('test_grid_R_%dn_%diter.csv', n, iter), result_matrix); 
#
#******************************************************************************
#{
elapsed_time = toc;
hours_time = floor(elapsed_time/3600);
minutes_time = floor((elapsed_time - (hours_time*3600))/60);
seconds_time = rem  ((elapsed_time - (hours_time*3600)),60);
elapsed_time_formated = [hours_time minutes_time seconds_time];
#}  
endfunction


% --------------------------------------------------------------
% Make one Seidel iteration
% --------------------------------------------------------------
function f = seidel(f, mask, na, mc)

# the area in the middle
for k = 2:mc
  for j = 2:na + (k-1)
    f(k,j) = f(k,j) + mask(k,j)*(0.25*( f(k-1,j) + f(k+1,j) + f(k,j-1) + f(k,j+1)) - f(k,j));
  end
end

% Symmetry on left boundary 

for k = 2:mc
  j = 1;
  f(k,j) = f(k,j) + mask(k,j)*(0.25*( f(k+1,j) + f(k-1,j) + f(k,j+1) + f(k,j+1)) - f(k,j));
end

% Symmetry on diagonal boundary   ??
for k = 2:mc
    j = na + (k-0);
  f(k,j) = f(k,j) + mask(k,j)*(0.25*( f(k+1,j) + f(k+1,j) + f(k,j-1) + f(k,j-1)) - f(k,j));
end
endfunction
%seidel = seidel(f, mask, na, mc)

% --------------------------------------------------------------
% Calculate E and J, budi mulyanto
% --------------------------------------------------------------
function I_1 = curr_dens_1(mc, na, f, delta_y, sigma, h) 
for k = 1:mc+0
  for j = 1:na+1+ (k-1) 
    E(k,j) = -(f(k+1,j) - f(k,j))/delta_y;
    J(k,j) = sigma*E(k,j);
  end
end

I_1 = 0;
for j = 1:na+1    
    I_upd = (J(1,j) * delta_y * h);
    I_1 = I_1 + I_upd;
endfor
endfunction

function I_2 = curr_dens_2(mc, na, f, delta_y, sigma, h) 
for k = 1:mc+0
  for j = 1:na+1+ (k-1) 
    E(k,j) = -(f(k+1,j) - f(k,j))/delta_y;
    J(k,j) = sigma*E(k,j);
  end
end

I_2 = 0;
  for j = 1:na+1+ (mc-1)    
    I_upd = (J(2,j) * delta_y*h);
    I_2 = I_2 + I_upd;
    end
endfunction

function I_mid = curr_dens_mid(mc, na, f, delta_y, sigma, h) 
for k = 1:mc+0
  for j = 1:na+1+ (k-1) 
    E(k,j) = -(f(k+1,j) - f(k,j))/delta_y;
    J(k,j) = sigma*E(k,j);
  end
end

I_mid = 0;
  for j = 1:na+1+ (round(mc/2)-1)    
    I_upd = (J(round(mc/2),j) * delta_y*h);
    I_mid = I_mid + I_upd;
    end
endfunction

function I_3 = curr_dens_3(mc, na, f, delta_y, sigma, h) 
for k = 1:mc+0
  for j = 1:na+1+ (k-1) 
    E(k,j) = -(f(k+1,j) - f(k,j))/delta_y;
    J(k,j) = sigma*E(k,j);
  end
end

I_old = 0;
I_avg = 0;
for k = 1:mc+0
  for j = 1:na+1+ (k-1)
    I_upd = (J(k,j) * delta_y*h);
    I_old = I_old + I_upd; # total I for each k
  end
    I_avg = (0 + I_old)/(k+0);
    k;
end
I_3 = I_avg/1;
endfunction

function E = E_field_1(mc, na, f, delta_y, sigma, n, iter) 
for k = 1:mc+0
  for j = 1:na+1+ (k-1) 
    E(k,j) = -(f(k+1,j) - f(k,j))/delta_y;
  end
end

#result_matrix = ([E]);
#dlmwrite (sprintf('test_grid_E_%dn_%diter.csv', n, iter), result_matrix);
endfunction

curr(40*3)















#{
function I_1 = current_1(curr_dens, na, mc, delta_y, sigma, h )
I_1 = 0;
for j = 1:na+1    
    I_upd = (curr_dens(1,j) * delta_y*h);
    I_1 = I_1 + I_upd;
endfor
endfunction


function I_2 = current_2(curr_dens, na, mc, delta_y, sigma, h )
I_2 = 0;
  for j = 1:na+1+ (mc-1)    
    I_upd = (curr_dens(mc,j) * delta_y*h);
    I_2 = I_2 + I_upd;
    end
endfunction
#}

#{
%function seidel_new = seidel_new(seidel, na, mc, delta_y, sigma, h)
function seidel_2 = seidel_2(seidel)
  for iter = 1:1000;
    f = seidel(f, mask, na, mc);
    n_iter = iter;
  endfor  
  %seidel_new = f;
endfunction
#}

#{
function f_iter = I_iteration(current_1, current_2, na, mc, delta_y, sigma, h)
I_1_old = 0;  
I_2_old = 0;  

for iter = 1:1000;
  f_iter = seidel(f, mask, na, mc);
  I_1 = current_1(curr_dens, na, mc, delta_y, sigma, h );
  I_2 = current_2(curr_dens, na, mc, delta_y, sigma, h );
  
  I_1_old = I_1;
    
  if (abs(I_1 - I_1_old)/I_1)<tol
    break
  else
    I_1_old = I_1;
  endif  
endfor

for iter = 1:1000;
  f = seidel(f, mask, na, mc);
  I_1 = current_1(curr_dens, na, mc, delta_y, sigma, h );
  I_2 = current_2(curr_dens, na, mc, delta_y, sigma, h );
  
  if (abs(I_2 - I_2_old)/I_2)<tol
    break
  else
    I_2_old = I_2;
  endif
endfor 
endfunction
#}
