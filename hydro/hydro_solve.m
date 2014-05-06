function [tOut uOut] = hydro_solve(params, uIn)

  time_steps = params.t_steps;


  global state_t_prev;
  global state_n;
  global prev_t;
  global fps; 
  state_t_prev = 0;
  state_n = 0;
  prev_t = 0;
  fps = 0;
 
  global xt
  xt = [];
  clear xt;
  global xt;
  xt = xtimer();

  % ode45 options
  
  if params.InitialStep ~= []
    opts = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol, 'OutputFcn', @(t,u,flag) hydro_plot(t, u, flag, params), 'InitialStep', params.InitialStep);
  else
    opts = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol, 'OutputFcn', @(t,u,flag) hydro_plot(t, u, flag, params));
  end

  %  2 -- 2nd order hydrodynamics, see mmatica file
  % -2 -- 2nd order hydro with terms O(ux,uy)^2 set to zero before inversion. the nr code is not as optimized as rel 2nd order and might actually run slower
  %  0 -- 0th order hydrodynamics
  switch params.hydroOrder
    case 2
      eom = @eom2;
    case -2
      eom = @eom2nr;
    case 0
      eom = @eom0;
  end

  hydro_begin=tic;
  % run the ode solver
  %[tOut, uOut] = ode45(@(t,u) eom(t,u,params), time_steps, uIn, opts);

  if length(time_steps) == 2
    sol = ode45(@(t,u) eom(t,u,params), time_steps, uIn, opts);
    tOut = sol
    uOut = 0;
  else
    [tOut uOut] = ode45(@(t,u) eom(t,u,params), time_steps, uIn, opts);
  end

  disp(sprintf('computed hydro simulation in %fsec', toc(hydro_begin)));

  xt.stats();
  return


end


