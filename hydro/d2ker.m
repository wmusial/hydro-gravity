
function ker = d2ker(order)
  switch (order)
    case 2
      ker = [1; -2; 1];
    case 4
      ker = [-1/12; 4/3; -5/2; 4/3; -1/12];
    case 6
      ker = [1/90; -3/20; 3/2; -49/18; 3/2; -3/20; 1/90];
    case 8
      ker = [-1/560; 8/315; -1/5; 8/5; -205/72; 8/5; -1/5; 8/315; -1/560];
  end
end
