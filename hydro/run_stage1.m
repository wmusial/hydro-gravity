function hydro_run_stage1()

  % runs a low-resolution hydro simulation in full solution mode so that we can extract the automatic time steps chosen by ode45

  jobname = 'turb-small';
  parts = 40; % every 50
  t_beg = 0;
  t_end = 20000; 
  Nx = 64;
  L = 1500; % ref 1
  Temp = 1; % ref 100
  hydroOrder = 2;

  AbsTol = 1e-6;
  RelTol = 1e-8;
  ic = 'turb';

  datestring = datestr(clock, 'yy-mm-dd-HHMM');

  tOut = [];

  InitialStep = [];
  fname = sprintf('%s-%s-s1', datestring, jobname);


  params_lo = hydro_params(Temp, Nx, L, [t_beg t_end], hydroOrder);
  params_lo.plotHydro = true;
  params_lo.RelTol = AbsTol;
  params_lo.AbsTol = RelTol;
  params_lo.hydroIC = ic;
  params_lo.InitialStep = InitialStep;

  % initializes derivative kernels, plotting 
  hydro_init(params_lo);

  for p=1:parts

    tbeg = (p-1).*t_end/parts;
    tend = (p).*t_end/parts;
    fname_part = sprintf('%s-%s-s1-T%d-L%d-Nx%d-tbeg%d-tend%d', datestring, jobname, round(Temp), round(L), Nx, round(tbeg), round(tend));
    fname_part

    t_steps = [tbeg tend];

    params_lo_part = params_lo.copy();
    params_lo_part.t_steps = t_steps;

    if p == 1
      % initial conditions, or replace with your own
      uIn = hydro_ic(params_lo_part);
    end

    % run the ode solver
    sol = hydro_solve(params_lo_part, uIn);

    % save results
    %save(fullfile(pwd, fname_part), 'jobname', 'params_lo', 'params_lo_part', 'sol', '-v7.3');

    % evaluate initial conditions for next part
    uIn = deval(sol, tend);

    % store computed times
    tOut = [tOut sol.x(1:end-1)];

    if length(sol.x) > 3
      InitialStep = sol.x(end-1) - sol.x(end-2);
    else
      InitialStep = [];
    end

    % clear memory
    sol = [];
    clear sol;

    save(fullfile(pwd, 'partial'), 'tOut');
  end


  tOut1 = [tOut tend];
  tOut1 = tOut1(:);
  save(fullfile(pwd, fname), 'jobname', 'params_lo',  'tOut1', '-v7.3');


end
