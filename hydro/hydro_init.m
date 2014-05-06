function hydro_init(params)

  % ****************************************************************
  %  output info
  % ****************************************************************

  % summary
  disp(sprintf('Grid Size   : (%d, %d)', params.Nx, params.Ny));
  disp(sprintf('Grid Length : %f', params.Lx));
  disp(sprintf('Cell Length : %f', params.Dx));
  %disp(sprintf('Time steps  : %d',  params.t));
  %disp(sprintf('Hydro Size  : %6.2fMB', 3*params.Nx.*params.Ny.*5.*params.Tsteps.*64/8/1024/1024));



  % ****************************************************************
  %  prepare plotting
  % ****************************************************************

  if (params.plotHydro)
    % prepare figures for plotting
    close all;
    global cm;

    % load color maps
    load('cm_turb', 'cm_turb');
    cm = cm_turb;

    global cmSize;
    cmSize = 512;

    % figures
    fig  = figure; 
    set(fig, 'WindowStyle', 'docked');
    set(fig, 'ColorMap', cm);
    set(gcf, 'renderer', 'zbuffer');

    refresh();
    pause(0.001);
  end
  

  % ****************************************************************
  %  compute spectral derivative kernels
  % ****************************************************************

  d_init(params.Nx, params.Dx);
end
