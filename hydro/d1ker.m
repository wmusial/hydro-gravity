function ker = d1ker(order)
  switch (order)
    case 2
      ker = [-1/2; 0; 1/2];
    case 4
      ker = [1/12; -2/3; 0; 2/3; -1/12];
    case 6
      ker = [-1/60; 3/20; -3/4; 0; 3/4; -3/20; 1/60];
    case 8
      ker = [1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280];
  end
  ker = flipud(ker);
end
