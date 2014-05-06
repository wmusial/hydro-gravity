function stop = hydro_writepng(tt, uu, flag, params)
  % writes png from uu
  % indended to be passed into ode45 as a plotting function

  % tell ode45 to keep going
  stop = 0;

  
  % exit if no plot on
  if (~params.plotHydro)
    return
  end

  if strcmp(flag, 'done')
    return
  end

  Nx = params.Nx;
  Ny = params.Ny;
  Nt = Nx.*Ny;

  % how many time steps did ode45 return?
  tsteps = size(uu,2);
  global state_t_prev;
  global state_n; 

  global prev_t;
  global fps;

  if prev_t ~= 0
    elapsed = toc(prev_t);
    fps_new = tsteps ./ elapsed;
    if fps == 0
      fps = fps_new;
    else
      r = 0.9;
      fps = r .* fps_new + (1-r) .* fps;
    end
  end
  prev_t = tic;


  t_curr = tt(tsteps);
  disp(sprintf('%5d  t = %.4f, Dt = %.4f, Nt = %d, %.1fFPS %.1fFPH', state_n, t_curr, t_curr - state_t_prev, tsteps, fps, fps.*60));

  state_t_prev = t_curr;
  state_n = state_n + 1;


  

  % collect results
  u = uu(:, tsteps);
  t = tt(:, tsteps);

  % massage the shapes
  T  = reshape(u(0*Nt+1 : 1*Nt), Nx, Ny);
  ux = reshape(u(1*Nt+1 : 2*Nt), Nx, Ny);
  uy = reshape(u(2*Nt+1 : 3*Nt), Nx, Ny);

  
  uxdy = dy(ux);
  uydx = dx(uy);

  % compute vorticity
  w = uxdy - uydx;

  global max_init_w;
  if max_init_w == []
    max_init_w = 0;
  end

  if max_init_w == 0
    max_init_w = max(abs(w(:)));
  end
    max_init_w = max(abs(w(:)));


  data = (w + max_init_w) ./ 2 ./ max_init_w;

  global cm cmSize
  data = ind2rgb(gray2ind(data, cmSize*2), cm);

  fname = sprintf('img/im-%04d.png', tt);
  imwrite(data, fname);


end
