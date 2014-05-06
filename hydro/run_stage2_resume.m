function hydro_run_stage2_resume(fname2)

  % runs a high-resolution hydro simulation at equally spaced discrete time points

  %*********************************************
  % load stage 2 params
  %*********************************************
  load(fname2, 'jobname', 'params_hi', 's2_parts');

  t_beg      = 0;
  t_end      = params_hi.t_steps(end);

  Nx = params_hi.Nx;
  Nt = Nx*Nx*5;

  % timestep
  dt = 5;

  datestring = datestr(clock, 'yy-mm-dd-HHMM');

  % initializes derivative kernels, plotting 
  hydro_init(params_hi);




  %*********************************************
  % load the previous part params
  %*********************************************
  load(fname2, 'params_hi_part');

  part_t_beg = params_hi_part.t_steps(1);
  part_t_end = params_hi_part.t_steps(end);
  p_previous = part_t_end .* s2_parts ./ t_end;
  
  % load the actual data
  load(fname2, 'uOut');

  uIn = uOut(end,:);
  uOut = [];
  clear uOut;


  %*********************************************
  % run simulations in parts
  %*********************************************

  for p=p_previous+1:s2_parts

    tbeg = (p-1).*t_end/s2_parts;
    tend = (p).*t_end/s2_parts;
    fname = sprintf('%s-%s-s2-T%d-L%d-Nx%d-tbeg%d-tend%d', datestring, jobname, round(Temp), round(L), Nx, round(tbeg), round(tend));
    fname

    t_steps = [tbeg : dt : tend];

    params_hi_part = params_hi.copy();
    params_hi_part.t_steps = t_steps;


    % run the ode solver
    [tOut, uOut] = hydro_solve(params_hi_part, uIn);

    % save last frame as input for next part
    uIn = uOut(end, :);

    save(fname, 'jobname', 'params_hi', 'params_hi_part', 's2_parts', 'tOut', 'uOut', '-v7.3');

    % clear mem
    uOut = [];
    clear uOut;
  end

end
