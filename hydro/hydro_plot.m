function stop = hydro_plot(tt, uu, flag, params)
  % plotting utility

  % tell ode45 to keep going
  stop = 0;

  global xt;
  xt.tic('plot');
  
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

  % legacy...?
  if (length(u(:)) > 4*Nt)
    P1 = reshape(u(3*Nt+1 : 4*Nt), Nx, Ny);
    P2 = reshape(u(4*Nt+1 : 5*Nt), Nx, Ny);
  else
    P1 = zeros(Nx, Ny);
    P2 = zeros(Nx, Ny);
  end

  
  %Tdx  = dx(T);
  %Tdy  = dy(T);
  %uxdx = dx(ux);
  uxdy = dy(ux);
  uydx = dx(uy);
  %uydy = dy(uy);

  % back of the envelope check if we're in the hydro regime
  if false
    T1 = sqrt( (Tdx.^2+Tdy.^2) )./T;

    u1 = sqrt(uxdx.^2 + uxdy.^2);
    u2 = sqrt(uydx.^2 + uydy.^2);

    disp(sprintf('hydro-check: T: (%.0e, %.0e), u: (%.0e, %.0e), u2: (%.0e, %.0e)', mean(mean(T1)), max(max(T1)), mean(mean(u1)), max(max(u1)), mean(mean(u2)), max(max(u2)) ));
  end

  % compute vorticity
  w = uxdy - uydx;


  figure(1);

  plotsx = 2;
  plotsy = 2;

  % custom imagesc
  function imagecsc(data)
    global cmSize;
    cmSize = 512;
    maxdata = max(max(abs(data)));
    data = data ./ maxdata;
    image(data.*cmSize + cmSize);
  end

  % temperature
  subplot(plotsy, plotsx, 1);
  imageh(params, uu, 1, sprintf('T @ t = %.1f', t));

  % w
  subplot(plotsy, plotsx, 2);
  imageh(params, w, 1, 'w');

  % ux
  subplot(plotsy, plotsx, 3);
  imageh(params, uu, 2, 'ux');

  % uy
  subplot(plotsy, plotsx, 4);
  imageh(params, uu, 3, 'uy');

  pause(0.001);
  drawnow();
  refresh();

  xt.toc('plot');
end
