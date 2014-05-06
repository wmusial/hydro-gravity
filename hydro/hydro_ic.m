function [uIn] = hydro_ic(params)

  Nx = params.Nx;
  Ny = params.Ny;


  % ****************************
  % initial conditions
  % ****************************
  T  = params.T .* ones(Nx, Ny);
  ux = zeros(Nx, Ny); % + 0.0001;
  uy = zeros(Nx, Ny); % + 0.0001;
  P1 = zeros(Nx, Ny);
  P2 = zeros(Nx, Ny);

  % ********************************************************
  %     V O R T I C E S
  % ********************************************************

  % initial conditions 
  if strcmp(params.hydroIC, 'vort')
    % specify vortices data
    vp =  [0.5 0.4 0.1 0.5; 0.5 0.6 0.1 +0.5] % [x, y, width, strength]
    %vp = [0.5 0.5 0.2 0.1];
    for i =1:size(vp,1)
      v = vp(i,:);
      vamp = v(4);
      vwidth = v(3);
      vx = v(1);
      vy = v(2);

      cx = vx.*params.Nx;
      cy = vy.*params.Ny;

      gauss_grid_width = vwidth .* params.Nx;
      sig2 = 2.*(gauss_grid_width).^2;

      for x = 1:params.Nx
        for y = 1:params.Nx
          for ddx = -1:1
            for ddy = -1:1
              %T(x, y) = T(x,y) + params.hydroPeakAmp .* exp(-((x-cx)^2 + (y-cy)^2)/sig2);
              ux(x,y) = ux(x,y) -vamp.*(y+ddy*params.Ny-cy).*exp(-((x+ddx*params.Nx-cx)^2 + (y+ddy*params.Ny-cy)^2)/sig2);
              uy(x,y) = uy(x,y) +vamp.*(x+ddx*params.Nx-cx).*exp(-((x+ddx*params.Nx-cx)^2 + (y+ddy*params.Ny-cy)^2)/sig2);
            end
          end
        end
      end  % end of fors
    end % one more for
  end % end of if

  % ********************************************************
  %     E X P    P E A K
  % ********************************************************

  if strcmp(params.hydroIC, 'peak')
    peak_amp = 0.1;
    peak_relwidth = 0.3;
    cx = 0.5 * params.Nx + 0.5;
    cy = 0.5 * params.Ny + 0.5;
    cx2 = 0.5 * params.Nx + 0.5;
    cy2 = 0.6 * params.Ny + 0.5;
    gauss_grid_width = peak_relwidth .* params.Nx;

    sig2 = 2.*(gauss_grid_width).^2;
    for x = 1:params.Nx
      for y = 1:params.Nx
        T(x, y) = T(x,y) + peak_amp .* exp(-((x-cx)^2 + (y-cy)^2)/sig2);
      end
    end
  end
  %ux = 0.99999999999999;

  [gx, gy] = ndgrid(params.Cx/params.Lx, params.Cy/params.Ly);
  gx = 2.*pi.*gx;
  gy = 2.*pi.*gy;


  % ********************************************************
  %     R A N D O M
  % ********************************************************

  rng(42);
  if strcmp(params.hydroIC, 'rand')
    rand_freq = 3;
    % smooth random gradients
    N = 1000;
    for s=1:N
      phasex = rand(1)*2*pi;
      phasey = rand(1)*2*pi;
      ampx = rand(1);
      ampy = rand(1);
      freqx  = round(rand(1)*rand_freq);
      freqy  = round(rand(1)*rand_freq);
      ux = ux + 1/N.*ampx.*sin(freqx.*gx + phasex).*sin(freqy.*gy + phasey);

      phasex = rand(1)*2*pi;
      phasey = rand(1)*2*pi;
      ampx = rand(1)*params.hydroPureWaveAmp;
      freqx  = round(rand(1)*rand_freq);
      freqy  = round(rand(1)*rand_freq);
      uy = uy + 1/N.*ampx.*sin(freqx.*gx + phasex).*sin(freqy.*gy + phasey); %.*sin(freqy*gy + phasey);
    end
    maxx = max(max(abs(ux)));
    maxy = max(max(abs(ux)));
    disp(sprintf('max ux: %f  max uy: %f', maxx, maxy));

    if true
      normalization = 0.1;
      ux = ux / max([maxx maxy]) .* normalization;
      uy = uy / max([maxx maxy]) .* normalization;
    end
  end


  % ********************************************************
  %     T U R B U L E N T
  % ********************************************************

  % L = 40
  % T = 250,500
    
  rng(42+1);
  if strcmp(params.hydroIC, 'turb')

    if true % allan
      % Lx = 750; % total size = 1500
      % Ly = 750; % total size = 1500
      % Nx = 305; 
      % Ny = 305;
      % Nt = Nx*Ny; 

      % To = 1;   % Temp = (3/4pi)T = 3/4pi
      Nw = 5;  % Nw full oscillations of uy in x
      Nk = 4;   % Number of fourier modes contributing to noise


      % t0 = 0;
      % tf = 3000;



      dux = 0*gx;
      duy = 0*gy;
      for kx = 1:Nk
          for ky = 1:Nk
              dux = dux + rand(1) * cos( kx.* (gx-pi) + ky.* (gy-pi) + 2*pi*rand(1));
              duy = duy + rand(1) * cos( kx.* (gx-pi) + ky.* (gy-pi) + 2*pi*rand(1));
          end
      end
      dux = dux / max(dux(:));
      duy = duy / max(duy(:));


      ux =                    0.01 * dux;
      uy = 0.1 * cos(Nw*gx) + 0.01 * duy;
      max(max(ux))
      max(max(uy))
    else % thesis



      % turbulent flow
      ux = ux + 0.1 * sin(5*gy);%.*(1+0.2.*cos(2*gx));

      for ix=1:5
        for iy=1:5
          uy = uy + 0.01 .* rand(1) .* sin(ix.*gx + iy.*gy + 2.*pi.*rand(1));
        end
      end
    end

      % old
      % turbulent flow
      % originally run with:
      %  L = 40
      %  T = 250/500
      %  t_end = 2500
      % 
      % equivalent to: 
      %  L = 10000
      %  T = 1/2
      %  t_end = 62500
      % 
      % P1/P2 scale as Temp.^3  !!!
      %
      % ux = ux + 0.1 * sin(2*gy);%.*(1+0.2.*cos(2*gx));
      % uy = uy + 0.005.*sin(2*gx+rand(1)*cos(3*gy));
      % uy = uy + 0.003.*sin(1*gx+rand(1)*cos(2*gy));
      % uy = uy + 0.0002.*sin(5*gx+rand(1)*cos(4*gx));
  end


  % ********************************************************
  %     W A V E
  % ********************************************************

  if strcmp(params.hydroIC, 'wave')
    wave_amp = 0.1;
    wave_freq = 5;
    ux = ux + wave_amp.*sin(wave_freq.*gy);
  end



  % ****************************
  %  rescale temperature
  % ****************************

  % T = T .* (4 * pi / 3);


  % repackage initial conditions into a vector
  Nt = Nx.*Ny;
  uIn = zeros(Nt.*5, 1);
  uIn(0*Nt+1 : 1*Nt) = T(:);
  uIn(1*Nt+1 : 2*Nt) = ux(:);
  uIn(2*Nt+1 : 3*Nt) = uy(:);
  uIn(3*Nt+1 : 4*Nt) = P1(:);
  uIn(4*Nt+1 : 5*Nt) = P2(:);

end
