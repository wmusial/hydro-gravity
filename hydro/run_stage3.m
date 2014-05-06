 function hydro_run_stage3(fname1, fnames2)

  %*********************************************
  % load stage 2 params
  %*********************************************
  load(fnames2{1}, 'jobname', 'params_hi', 's2_parts');

  t_beg      = 0;
  t_end      = params_hi.t_steps(end);
  L          = params_hi.Lx;
  Temp       = params_hi.T;
  hydroOrder = params_hi.hydroOrder

  AbsTol = params_hi.AbsTol;
  RelTol = params_hi.RelTol;
  ic     = params_hi.hydroIC;

  Nx = params_hi.Nx;
  Nt = Nx*Nx*5;

  datestring = datestr(clock, 'yy-mm-dd-HHMM');

  % initializes derivative kernels, plotting 
  hydro_init(params_hi);








  %*********************************************
  % load stage1 timestep information
  %*********************************************
  load(fname1, 'tOut1', 'jobname');

  % compute time steps ode45 converged on while solving the pre-run
  tOut1 = tOut1(:);
  dtOut1 = conv(tOut1, [1, -1], 'valid');
  dtOut1 = medfilt2(dtOut1, [7; 1]);

  %*********************************************
  % metric evaluation points
  %*********************************************

  % metric evaluation point time -- subdivide the interval
  Nmep = 30;
  mep_time = [1:Nmep]/(Nmep) .* t_end;

  % ... or choose manually 
  % mep_time = 2000;
  % Nmep = length(mep_time);

  % time differencing order
  tord = 4;

  % initialize desired timesteps with zero
  mep_tsteps = zeros(Nmep, tord+1);
  mep_dt     = zeros(1, Nmep);
  for mi = 1:Nmep
    ti = find(tOut1<mep_time(mi), 1, 'last');
    mep_dt(mi) = dtOut1(ti);
    mep_tsteps(mi,:) = [-tord/2:tord/2].*mep_dt(mi)+mep_time(mi);
  end

  mep_time
  mep_dt

  %*********************************************
  % load stage2 params
  %*********************************************
  load(fnames2{1}, 'params_hi', 's2_parts'); % load part-independent params for stage 2
  Nx = params_hi.Nx;
  Nt = Nx*Nx*5;
  % simulation time domain
  t_beg = 0;
  t_end      = params_hi.t_steps(end);
  assert(length(fnames2) == s2_parts)

  %*********************************************
  % initialize u vars
  %*********************************************

  uOut3     = zeros(Nmep, tord+1, Nt);
  udtOut3   = zeros(Nmep, tord+1, Nt);
  udtdtOut3 = zeros(Nmep, tord+1, Nt);
  Delta = zeros(Nmep, tord+1, Nx, 5*Nx);


  datestring = datestr(clock, 'yy-mm-dd-HHMM');


  %*********************************************
  % Simulations
  %*********************************************

  hydro_init(params_hi);

  for fi = 1:length(fnames2)
    fprintf('file: %s\n', fnames2{fi});
    load(fnames2{fi}, 'params_hi_part') % load part-dependent params_hi_part
    fprintf('time range: (%.1f, %.1f)\n', params_hi_part.t_steps(1), params_hi_part.t_steps(end));


    mep_toevaluate = [];
    % run simulations for metric points within current file timespan
    for mi=1:Nmep

      % skip mep if outside of file's time interval
      min_t = min(mep_tsteps(mi, :));
      if (min_t <= params_hi_part.t_steps(1)) | (min_t > params_hi_part.t_steps(end))
        continue;
      end

      % flag mep for evaluation
      mep_toevaluate = [mep_toevaluate mi];
    end

    % skip the current file if no metric points 
    if length(mep_toevaluate) == 0
      disp('no metric points to evaluate, skipping\n');
      continue;
    end

    % if evaluating metric points, load the file
    load(fnames2{fi}, 'tOut2', 'uOut2');

    for mi=mep_toevaluate
      ti = find(tOut2 < min(mep_tsteps(mi, :)), 1, 'last');
      uIn = uOut2(ti, :);

      tsteps = [tOut2(ti), mep_tsteps(mi, :)];

      params_hi_part.t_steps = tsteps;

      [tMep, uMep] = hydro_solve(params_hi_part, uIn);
      % delete the first t=0 point 
      uMep(1,:) = [];
      tMep(1)   = [];

      % initialize the derivative matrices
      udtMep   = 0 * uMep;
      udtdtMep = 0 * uMep;

      
      

      % evaluate derivatives
      for ti=1:length(tMep)
        u = uMep(ti,:);
        udt   = eom2 (tMep(ti), u, params_hi_part);
        udtdt = deom2(tMep(ti), u, udt, params_hi_part);
        %size([a b c d e]);
        udtOut3  (mi, ti, :) = udt;
        udtdtOut3(mi, ti, :) = udtdt;
        if true:
          [a b c d e] = constituent_relations(tMep(ti), u, udt, udtdt, params_hi_part);
          Delta(mi, ti, :,:) = [a b c d e];
        else:
          Delta = []
        end
      end
      uOut3(mi, :, :) = uMep;
    end
  end

  uOut3    (:,:,Nx*Nx*3+1:Nx*Nx*5) = [];
  udtOut3  (:,:,Nx*Nx*3+1:Nx*Nx*5) = [];
  udtdtOut3(:,:,Nx*Nx*3+1:Nx*Nx*5) = [];


  fname = sprintf('%s-%s-s3-mep%d', datestring, jobname, Nmep);
  save(fname, 'jobname', 'params_hi', 'uOut3', 'udtOut3', 'udtdtOut3', 'Delta', 'mep_tsteps', 'mep_time', 'mep_dt', '-v7.3');
end
