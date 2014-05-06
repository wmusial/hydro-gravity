classdef hydro_params < handle

  properties (Constant)
    A=sqrt(3)*pi/18 - log(3)/2;
  end

  % parameters that are cast in stone are initialized with the constructor
  properties (GetAccess='public', SetAccess='private')

    % background T offset
    T = 1

    % number of samples in each dimension
    Nx        = 32
    Ny        = 32

    % physical size in each dimension
    Lx        = 100
    Ly        = 100

    % cell coordinate values
    Cx
    Cy

    % coordinate distance between adj cells
    Dx
    Dy

    % order of hydro equations
    hydroOrder = 2;

  end

  properties

    % time steps 
    t_steps = 0;

    plotHydro = false;


    hydroIC = 'turb';

    % tolerances
    RelTol = 1e-8;
    AbsTol = 1e-10;
    InitialStep = []; % only used for the full solution-mode ode45 run
    % turb | rand | wave | peak | vort
    %
    % turb - turbulent flow
    % rand - random 
    % wave - pure sin wave
    % peak - gaussian peak
    % vort - vortices
    % 
    % see hydro_ic.m for details
   
  end
  methods 
    %function this=hydro_params()
    function this=hydro_params(varargin)
      % T, Nx, L, t_end, hydroorder, timeorder


      this.T = varargin{1};
      
      % number of samples
      this.Nx = varargin{2};
      this.Ny = varargin{2};

      % physical length
      this.Lx = varargin{3};
      this.Ly = varargin{3};

      % time limits
      %this.t_beg = 0;
      %this.t_end = varargin{4};

      this.t_steps = varargin{4};

      this.hydroOrder = varargin{5};

      % cell values
      this.Cx = (0:this.Nx-1)/this.Nx*this.Lx;
      this.Cy = (0:this.Ny-1)/this.Ny*this.Ly;

      % separation between cells
      this.Dx = this.Lx / this.Nx;
      this.Dy = this.Ly / this.Ny;

    end
    function r=copy(this)
       % provide a "copy constructor" to circumvent matlab's OO
      r = hydro_params(this.T, this.Nx, this.Lx, this.t_steps, this.hydroOrder);
      r.plotHydro = this.plotHydro;
      r.hydroIC = this.hydroIC;
      r.RelTol = this.RelTol;
      r.AbsTol = this.AbsTol;
      r.InitialStep = this.InitialStep;
    end
  end % end of methods
end
