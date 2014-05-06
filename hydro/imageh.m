% custom imagesc for hydro
function imageh(params, u, n, varargin)
  Nx = params.Nx;
  Ny = params.Ny;
  Nt = Nx.*Ny;

  if n == 0 % vorticity
    n = 2;
    ux = reshape(u((n-1)*Nt+1:n*Nt), Nx, Ny);
    n = 3;
    uy = reshape(u((n-1)*Nt+1:n*Nt), Nx, Ny);
    
    data = dy(ux) - dx(uy);
    n = 0;
  else
    data = reshape(u((n-1)*Nt+1:n*Nt), Nx, Ny);
  end
  %data = log(abs(dy(dx(data))));

    %data = data - mean(data(:));
  offset = false;
  if mean(data(:)) > 0.5 .* std(data(:))
    data = data - mean(data(:));
    offset = true;
  end

  % build title string
  tstring = '';
  if (length(varargin) >= 1) & (~strcmp(varargin{1}, ''))
    tstring = [varargin{1} '   ---   '];
  end

  % take varargin{2} dx derivatives
  if length(varargin) >= 2
    for i=1:varargin{2}
      data = dx(data);
    end
  end
  % take varargin{2} dy derivatives
  if length(varargin) >= 3
    for i=1:varargin{2}
      data = dy(data);
    end
  end

  % build title string
  tstring = [tstring  sprintf('max: %.4g, avg: %.4g, std: %.4g', max(abs(data(:))), mean(data(:)), std(data(:)))];

  if offset
    tstring = [tstring ', offset'];
  end

  global cmSize;
  maxdata = max(max(abs(data)));
  data = data ./ maxdata;
  image(data.*cmSize + cmSize);
  axis square;


  title(tstring);


end
