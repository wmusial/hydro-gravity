function hplayer(filename)

  load(filename, 'params_hi', 'tOut', 'uOut');
  params  = params_hi
  params.plotHydro = true;

  % initialize derivative kernels
  hydro_init(params);

  Nx = params.Nx;
  Ny = params.Ny;
  Nt = Nx.*Ny

  for t=1:length(tOut(:))
    hydro_writepng(t+1809, reshape(uOut(t,:), Nt*5, 1),'', params);
    pause(0.01);
  end
end


















