function hydro_run_stage2(fname1)

  % runs a high-resolution hydro simulation at equally spaced discrete time points

  Nx = 306;
  Nt = Nx*Nx*5;
  s2_parts   = 3;


  %*********************************************
  % load stage1 timestep information
  %*********************************************
  load(fname1, 'params_lo', 'jobname');

  t_beg      = 0;
  t_end      = params_lo.t_steps(end);
  L          = params_lo.Lx;
  Temp       = params_lo.T;
  hydroOrder = params_lo.hydroOrder

  AbsTol = params_lo.AbsTol;
  RelTol = params_lo.RelTol;
  ic     = params_lo.hydroIC;

  % timestep
  dt = 100;

  datestring = datestr(clock, 'yy-mm-dd-HHMM');




  %*********************************************
  % initialize hi-res params
  %*********************************************

  t_steps = [t_beg, t_end]; % this is to be overriden in the loop
  params_hi = hydro_params(Temp, Nx, L, t_steps, hydroOrder);
  params_hi.plotHydro = true;
  params_hi.RelTol = RelTol;
  params_hi.AbsTol = AbsTol;
  params_hi.hydroIC = ic;

  % initializes derivative kernels, plotting 
  hydro_init(params_hi);

  % initial conditions
  uIn = hydro_ic(params_hi);

  %*********************************************
  % run simulations in parts
  %*********************************************

  for p=1:s2_parts

    tbeg = (p-1).*t_end/s2_parts;
    tend = (p).*t_end/s2_parts;
    fname = sprintf('%s-%s-s2-T%d-L%d-Nx%d-tbeg%d-tend%d', datestring, jobname, round(Temp), round(L), Nx, round(tbeg), round(tend));
    fname

    t_steps = [tbeg : dt : tend];

    params_hi_part = params_hi.copy();
    params_hi_part.t_steps = t_steps;


    % run the ode solver
    [tOut2, uOut2] = hydro_solve(params_hi_part, uIn);

    save(fname, 'jobname', 'params_hi', 'params_hi_part', 's2_parts', 'tOut2', 'uOut2', '-v7.3');

    % keep the last frame as input for next part
    uIn = uOut2(end, :);

    % clear mem
    uOut2 = [];
    clear uOut2;
  end


end
