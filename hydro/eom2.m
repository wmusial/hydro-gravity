function udt = eom2(t, u, params)
  drawnow();
%   choice of Pi: 
%   Pi_xx =  2 P1 (1+UX^2) +  2 UX UY P2
%   Pi_yy = -2 P1 (1+UY^2) +  2 UX UY P2
%   Pi_xy =  P2 ( 2 + UX^2 + UY^2 )

    A = params.A;
    A1 = A+1;
    A2 = 2.*A;
    Nx = params.Nx;
    Ny = params.Ny;
    Nt = Nx.*Ny;

    %u = gpuArray(uu);
    


    % unpack variables 
    T  = reshape(u(0*Nt+1 : 1*Nt), Nx, Ny);
    ux = reshape(u(1*Nt+1 : 2*Nt), Nx, Ny);
    uy = reshape(u(2*Nt+1 : 3*Nt), Nx, Ny);
    P1 = reshape(u(3*Nt+1 : 4*Nt), Nx, Ny);
    P2 = reshape(u(4*Nt+1 : 5*Nt), Nx, Ny);

    ut = sqrt(1+ux.^2+uy.^2);
    
    global xt;

    xt.tic('eom2 derivs');

    Tdx  = dx(T);
    Tdy  = dy(T);
    uxdx = dx(ux);
    uxdy = dy(ux);
    uydx = dx(uy);
    uydy = dy(uy);
    P1dx = dx(P1);
    P1dy = dy(P1);
    P2dx = dx(P2);
    P2dy = dy(P2);

    xt.toc('eom2 derivs');
    xt.tic('eom2 precomp');

    utdx = (1/2).*ut.^(-1).*(2.*ux.*uxdx+2.*uy.*uydx);
    utdy = (1/2).*ut.^(-1).*(2.*ux.*uxdy+2.*uy.*uydy);


    ut2 = ut.^2;
    ut3 = ut.^3;
    ut4 = ut.^4;
    ux2 = ux.^2;
    ux3 = ux.^3;
    ux4 = ux.^4;
    uy2 = uy.^2;
    uy3 = uy.^3;
    uy4 = uy.^4;
    T2 = T.^2;
    T3 = T.^3;
    T4 = T.^4;

    xt.toc('eom2 precomp');

    xt.tic('eom2 M');


    m11 = (-3).*T2.*ut2;
    m12 = (-3/2).*T3.*ux+P2.*uy.*((-2)+ux2+(-1).*uy2)+(-2).*P1.*ux.*(1+uy2);
    m13 = (1/2).*((-3).*T3+4.*P1.*(1+ux2)).*uy+P2.*ux.*((-2)+(-1).*ux2+uy2);

    m22 = (1/4).*ut.^(-2).*ux.*uy.*(3.*T3.*ut2.*(2+ux2+2.*uy2)+2.*P2.*ux.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+4.*P1.*(ux4.*(1+(-3).*uy2)+2.*(1+uy2)+(-3).*ux2.*((-1)+uy4)));
    m23 = (1/4).*ut.^(-2).*ux.*uy.*(ux.*uy.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2))+2.*P2.*(4+(-3).*ux4+(-3).*ux.^6+6.*uy2+4.*uy4+ux2.*(4+3.*(uy2+uy4))));
    m24 = 2.*ux2.*(1+ux2).*uy;
    m25 = ux.*uy2.*(2+3.*ux2+uy2);

    m32 = (1/4).*ut.^(-2).*ux.*uy.*((-1).*ux.*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2))+2.*P2.*(3.*ux2.*(2+uy2)+ux4.*(4+3.*uy2)+(-1).*(1+uy2).*((-4)+3.*uy4)));
    m33 = (1/4).*ut.^(-2).*ux.*uy.*(3.*T3.*ut2.*(2+2.*ux2+uy2)+4.*P1.*((-2).*(1+ux2)+3.*((-1)+ux4).*uy2+((-1)+3.*ux2).*uy4)+2.*P2.*ux.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4)));
    m34 = (-2).*ux.*uy2.*(1+uy2);
    m35 = ux2.*uy.*(2+ux2+3.*uy2);

    m42 = (1/12).*(3.*T3.*(T3.*ux.*(2+ux2+3.*uy2)+A2.*uy.*(2.*P1.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-1).*P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)))+2.*A1.*(8.*P1.^2.*ux.*(1+uy2).*(2+ux2+uy2)+6.*P2.*T3.*uy.*(2+ux4+3.*uy2+uy4+ux2.*(1+(-2).*uy2))+P1.*(3.*T3.*ux.*(6+7.*uy2+4.*uy4+ux2.*(3+(-4).*uy2))+4.*P2.*uy.*((-1).*ux4+(2+uy2).^2))));
    m43 = (1/12).*(3.*T3.*((-1).*T3.*uy.*(2+3.*ux2+uy2)+A2.*ux.*(2.*P1.*ux.*uy.*((-1).*ux2+uy2)+P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)))+2.*A1.*((-8).*P1.^2.*(1+ux2).*uy.*(2+ux2+uy2)+(-6).*P2.*T3.*ux.*(2+ux4+uy2+uy4+ux2.*(3+(-2).*uy2))+P1.*(4.*P2.*ux.*((2+ux2).^2+(-1).*uy4)+3.*T3.*uy.*(ux2.*(7+4.*ux2+(-4).*uy2)+3.*(2+uy2)))));
    m44 = A1.*T3.*ut2.*(2+ux2+uy2);

    m52 = (1/6).*ux.*uy.*(3.*T3.*uy.*(T3.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)))+A1.*(4.*P1.*(1+uy2).*((-3).*T3.*(1+2.*ux2).*uy+2.*P2.*ux.*(2+ux2+uy2))+P2.*(4.*P2.*uy.*((-1).*ux4+(2+uy2).^2)+3.*T3.*ux.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2)))));
    m53 = (1/6).*ux.*uy.*(3.*T3.*ux.*(T3.*(1+ux2)+A2.*(P2.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1.*(1+ux2).*(1+uy2)))+A1.*(12.*P1.*T3.*ux.*(1+ux2).*(1+2.*uy2)+4.*P2.^2.*ux.*((2+ux2).^2+(-1).*uy4)+P2.*uy.*((-8).*P1.*(1+ux2).*(2+ux2+uy2)+3.*T3.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2))));
    m55 = A1.*T3.*ut2.*ux.*uy.*(2+ux2+uy2);


    xt.toc('eom2 M');
    xt.tic('eom2 R');



    R1 = (1/2).*((-6).*T2.*ut.*(Tdx.*ux+Tdy.*uy)+(-3).*T3.*ut.*(uxdx+ut.*(utdx.*ux+utdy.*uy)+(-1).*ux.*(ux.*uxdx+uy.*(uxdy+uydx))+uydy+(-1).*uy2.*uydy)+2.*(2.*P1.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy))+P2.*(utdx.*uy.*(2+uy2)+(-1).*ut.*(1+ut2).*(uxdy+uydx)+ux.*(utdy.*(2+ux2+3.*uy2)+uy.*(3.*utdx.*ux+(-2).*ut.*(uxdx+uydy))))));
    R2 = (1/4).*((-8).*P1.*uy2.*(utdy+utdy.*uy2+(-1).*ut.*uy.*uydy)+ux.*(2.*uy.*(4.*P1dx.*ut+3.*T2.*Tdx.*ut+2.*P2dy.*(ut+ut3)+8.*P2.*utdy+uy.*((4.*P1+3.*T3).*ut.*uxdy+4.*P2.*(2.*utdy.*uy+ut.*(uxdx+(-1).*uydy))))+ux.*(8.*P1.*utdy+uy.*(2.*(4.*P2dx.*ut.*uy+P2.*(3.*utdx.*uy.*(2+uy2)+ut.*((-10).*uxdy+(-2).*uydx+(-3).*uy2.*(uxdy+uydx))))+(-3).*T3.*ut.*((-1).*uxdx+3.*ut.*utdy.*uy+uydy+(-3).*uy2.*uydy)+4.*P1.*((-3).*utdy.*uy.*(1+uy2)+ut.*(uxdx+uydy+3.*uy2.*uydy)))+ux.*(P1.*((-8).*ut.*uxdy+12.*utdx.*uy+ux.*(8.*utdy+12.*utdx.*ux.*uy+(-12).*ut.*uxdx.*uy))+uy.*(8.*P1dx.*ut+(-9).*T3.*ut.*(ut.*utdx+(-1).*ux.*uxdx+(-1).*uy.*(uxdy+uydx))+2.*P2.*(utdy.*(14+3.*ux2+9.*uy2)+3.*(ux.*(3.*utdx.*uy+(-1).*ut.*(uxdy+uydx))+(-2).*ut.*uy.*(uxdx+uydy))))))));
    R3 = (1/4).*((-8).*P1.*uy2.*(utdx+utdx.*uy2+(-1).*ut.*uy.*uydx)+ux.*(uy.*(6.*T2.*Tdy.*ut+4.*ut.*((-2).*P1dy+P2dx+P2dx.*ut2)+16.*P2.*utdx+uy.*((-2).*uy.*(4.*P1dy.*ut+P2.*(3.*ut.*uxdy.*uy+utdx.*((-14)+(-3).*uy2)))+(-2).*P2.*ut.*(2.*uxdy+(7+3.*ut2).*uydx)+(-3).*T3.*ut.*(uxdx+3.*ut.*utdy.*uy+((-1)+(-3).*uy2).*uydy)+(-4).*P1.*(3.*utdy.*uy.*(1+uy2)+ut.*(uxdx+uydy+(-3).*uy2.*uydy))))+ux.*(4.*P1.*(utdx.*(1+ux2).*(2+3.*uy2)+(-1).*ut.*(ux.*uxdx.*(2+3.*uy2)+2.*uy.*uydx))+uy.*((-1).*ut.*((-8).*P2dy.*uy+3.*T3.*((-3).*ux.*uxdx.*uy+(-2).*uydx+(-3).*uy.*((-1).*ut.*utdx+uy.*(uxdy+uydx))))+2.*P2.*(3.*utdy.*uy.*(2+ux2+3.*uy2)+utdx.*ux.*(8+9.*uy2)+ut.*((-4).*uxdx+(-3).*ux.*uxdy.*uy+4.*uydy+(-6).*uy2.*(uxdx+uydy)))))));
    R4 = (1/12).*(3.*A2.*T3.*((-2).*P1.*ux.*uy.*(ux2+(-1).*uy2)+P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2)).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+6.*T4.*ut.*(2.*P1.*(1+ut2)+T2.*(uxdx+(-1).*uydy))+3.*T.^6.*uy.*(utdy.*(2+3.*uy2+uy4)+ut.*uy.*(uxdx+(-1).*(3+uy2).*uydy))+(-3).*T.^6.*ux.*(utdx.*(2+3.*ux2+ux4+uy2+(-1).*uy4)+ut.*uy3.*(uxdy+uydx)+ux.*(utdy.*((-1)+ux2).*uy+(-1).*ut.*(uxdx.*(3+ux2+(-1).*uy2)+ux.*uy.*(uxdy+uydx)+((-1)+uy2).*uydy)))+2.*A1.*(6.*P1dx.*T3.*ut.*(1+ut2).*ux+(-8).*P1.^2.*(1+ut2).*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy))+6.*T3.*(ut.*uy.*(P1dy.*(1+ut2)+P2.*uxdy.*uy.*(2+uy2))+P2.*ux.*(uy.*((-2).*utdy.*uy3+ut.*(uxdx.*(2+uy2)+((-2)+uy2).*uydy))+ux.*(2.*utdy.*ux.*uy2+(-2).*utdx.*uy.*((-1).*ux2+uy2)+(-1).*ut.*(uy2.*(uxdy+(-1).*uydx)+2.*uydx+ux.*(ux.*uydx+uy.*(uxdx+uydy))))))+P1.*(4.*P2.*(1+ut2).*((-1).*utdx.*uy.*(2+3.*ux2+uy2)+(-1).*utdy.*ux.*(2+ux2+3.*uy2)+ut3.*(uxdy+uydx)+ut.*(uxdy+uydx+2.*ux.*uy.*(uxdx+uydy)))+T3.*(utdx.*ux.*(12+ux2.*(30+18.*ux2+12.*uy2)+uy2.*(30+18.*uy2))+3.*(2.*utdy.*uy.*(2+3.*ux4+5.*uy2+3.*uy4+ux2.*(5+2.*uy2))+(-1).*ut.*((-6).*uxdx+(-6).*uydy+uy2.*((-3).*uxdx+uydy+6.*uy2.*uydy)+ux.*(2.*uy.*(uxdy.*(2+uy2)+(2+3.*uy2).*uydx)+ux.*(uxdx+2.*ux.*(3.*ux.*uxdx+uy.*(3.*uxdy+uydx))+(-3).*uydy+2.*uy2.*(uxdx+uydy)))))))));
    R5 = (1/6).*(3.*A2.*T3.*ux.*uy.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)).*((-1).*utdy.*ux+utdx.*uy+ut.*(uxdy+(-1).*uydx))+3.*T4.*ut.*ux.*uy.*(2.*P2.*(1+ut2)+T2.*(uxdy+(-1).*ut.*(utdx.*(1+ux2).*uy+utdy.*ux.*(1+uy2))+uydx+uy2.*(uxdy+uydx)+ux.*(ux.*(uxdy+ux.*uxdx.*uy+uydx+uy2.*(uxdy+uydx))+uy.*(uxdx+uydy+uy2.*uydy))))+A1.*((-12).*P1.*T3.*uy3.*(utdy+utdy.*uy2+(-1).*ut.*uy.*uydy)+ux.*(uy.*(6.*P2dy.*T3.*ut.*(1+ut2).*uy+4.*P2.^2.*(1+ut2).*((-1).*utdx.*uy.*(2+uy2)+ut.*(1+ut2).*(uxdy+uydx))+3.*P2.*T3.*(6.*ut.*(uxdx+uydy)+uy.*(utdy.*(12+18.*uy2+4.*uy4)+ut.*uy.*(3.*uxdx+((-9)+(-4).*uy2).*uydy)))+4.*P1.*((-3).*T3.*uy.*(utdx+utdx.*uy2+(-1).*ut.*uy.*uydx)+2.*P2.*(1+ut2).*(utdy.*uy.*(1+uy2)+ut.*(uxdx+(-1).*(1+uy2).*uydy))))+ux.*(4.*P1.*((-2).*P2.*(1+ut2).*(utdx+utdx.*ux2+(-1).*ut.*ux.*uxdx).*uy+3.*T3.*(ux.*(utdx+ux.*(utdx.*ux+(-1).*ut.*uxdx+utdy.*uy).*(1+uy2)+uy.*((-1).*utdx.*uy3+(-1).*ut.*(uxdy+uy2.*(uxdy+(-1).*uydx))))+uy.*(utdy+(-1).*utdy.*uy4+ut.*uy3.*uydy)))+uy.*(6.*P2dx.*T3.*ut.*(1+ut2)+P2.*((-4).*P2.*(1+ut2).*(utdy.*(2+ux2+3.*uy2)+uy.*(3.*utdx.*ux+(-2).*ut.*(uxdx+uydy)))+T3.*(utdx.*(36+uy2.*(54+12.*uy2)+ux2.*(54+12.*ux2+48.*uy2))+(-12).*ut.*uy.*(uxdy.*(3+2.*uy2)+(3+uy2).*uydx)+3.*ux.*(2.*utdy.*uy.*(9+2.*ux2+8.*uy2)+ut.*((-9).*uxdx+(-4).*ux.*(ux.*uxdx+uy.*(uxdy+2.*uydx))+3.*uydy+(-8).*uy2.*(uxdx+uydy))))))))));


    xt.toc('eom2 R');
    xt.tic('eom2 mdet');


    mdet = m11.*m25.*m34.*m43.*m52+(-1).*m11.*m24.*m35.*m43.*m52+(-1).*m11.*m25.*m33.*m44.*m52+m11.*m23.*m35.*m44.*m52+(-1).*m11.*m25.*m34.*m42.*m53+m11.*m24.*m35.*m42.*m53+m11.*m25.*m32.*m44.*m53+(-1).*m11.*m22.*m35.*m44.*m53+(-1).*m11.*m24.*m33.*m42.*m55+m11.*m23.*m34.*m42.*m55+m11.*m24.*m32.*m43.*m55+(-1).*m11.*m22.*m34.*m43.*m55+(-1).*m11.*m23.*m32.*m44.*m55+m11.*m22.*m33.*m44.*m55;

   
    i11 = m35.*((-1).*m24.*m43.*m52+m23.*m44.*m52+m24.*m42.*m53+(-1).*m22.*m44.*m53)+m25.*(m34.*m43.*m52+(-1).*m33.*m44.*m52+(-1).*m34.*m42.*m53+m32.*m44.*m53)+((-1).*m24.*m33.*m42+m23.*m34.*m42+m24.*m32.*m43+(-1).*m22.*m34.*m43+(-1).*m23.*m32.*m44+m22.*m33.*m44).*m55;
    i12 = (-1).*m13.*(m35.*m44.*m52+m34.*m42.*m55+(-1).*m32.*m44.*m55)+m12.*(m35.*m44.*m53+m34.*m43.*m55+(-1).*m33.*m44.*m55);
    i13 = m25.*m44.*(m13.*m52+(-1).*m12.*m53)+(m13.*m24.*m42+(-1).*m12.*m24.*m43+(-1).*m13.*m22.*m44+m12.*m23.*m44).*m55;
    i14 = (-1).*(m25.*m34+(-1).*m24.*m35).*(m13.*m52+(-1).*m12.*m53)+((-1).*m13.*m24.*m32+m12.*m24.*m33+m13.*m22.*m34+(-1).*m12.*m23.*m34).*m55;
    i15 = (m25.*m34+(-1).*m24.*m35).*(m13.*m42+(-1).*m12.*m43)+((-1).*m13.*m25.*m32+m12.*m25.*m33+m13.*m22.*m35+(-1).*m12.*m23.*m35).*m44;
    
    i22 = (-1).*m11.*(m35.*m44.*m53+m34.*m43.*m55+(-1).*m33.*m44.*m55);
    i23 = m11.*(m25.*m44.*m53+m24.*m43.*m55+(-1).*m23.*m44.*m55);
    i24 = m11.*((-1).*m25.*m34.*m53+m24.*m35.*m53+(-1).*m24.*m33.*m55+m23.*m34.*m55);
    i25 = m11.*(m25.*m34.*m43+(-1).*m24.*m35.*m43+(-1).*m25.*m33.*m44+m23.*m35.*m44);
    
    i32 = m11.*(m35.*m44.*m52+m34.*m42.*m55+(-1).*m32.*m44.*m55);
    i33 = (-1).*m11.*(m25.*m44.*m52+m24.*m42.*m55+(-1).*m22.*m44.*m55);
    i34 = m11.*(m25.*m34.*m52+(-1).*m24.*m35.*m52+m24.*m32.*m55+(-1).*m22.*m34.*m55);
    i35 = m11.*((-1).*m25.*m34.*m42+m24.*m35.*m42+m25.*m32.*m44+(-1).*m22.*m35.*m44);

    i42 = m11.*((-1).*m35.*m43.*m52+m35.*m42.*m53+(-1).*m33.*m42.*m55+m32.*m43.*m55);
    i43 = m11.*(m25.*m43.*m52+(-1).*m25.*m42.*m53+m23.*m42.*m55+(-1).*m22.*m43.*m55);
    i44 = m11.*((-1).*m25.*m33.*m52+m23.*m35.*m52+m25.*m32.*m53+(-1).*m22.*m35.*m53+(-1).*m23.*m32.*m55+m22.*m33.*m55);
    i45 = m11.*(m25.*m33.*m42+(-1).*m23.*m35.*m42+(-1).*m25.*m32.*m43+m22.*m35.*m43);

    i52 = m11.*(m34.*m43.*m52+(-1).*m33.*m44.*m52+(-1).*m34.*m42.*m53+m32.*m44.*m53);
    i53 = m11.*((-1).*m24.*m43.*m52+m23.*m44.*m52+m24.*m42.*m53+(-1).*m22.*m44.*m53);
    i54 = m11.*(m24.*m33.*m52+(-1).*m23.*m34.*m52+(-1).*m24.*m32.*m53+m22.*m34.*m53);
    i55 = m11.*((-1).*m24.*m33.*m42+m23.*m34.*m42+m24.*m32.*m43+(-1).*m22.*m34.*m43+(-1).*m23.*m32.*m44+m22.*m33.*m44);

    Tdt  = -(i11.*R1 + i12.*R2 + i13.*R3 + i14.*R4 + i15.*R5) ./ mdet;
    uxdt = -(          i22.*R2 + i23.*R3 + i24.*R4 + i25.*R5) ./ mdet;
    uydt = -(          i32.*R2 + i33.*R3 + i34.*R4 + i35.*R5) ./ mdet;
    P1dt = -(          i42.*R2 + i43.*R3 + i44.*R4 + i45.*R5) ./ mdet;
    P2dt = -(          i52.*R2 + i53.*R3 + i54.*R4 + i55.*R5) ./ mdet;


    xt.toc('eom2 mdet');


    udt = zeros(Nx.*Ny.*5, 1);

    udt(0*Nt+1 : 1*Nt) = Tdt(:);
    udt(1*Nt+1 : 2*Nt) = uxdt(:);
    udt(2*Nt+1 : 3*Nt) = uydt(:);
    udt(3*Nt+1 : 4*Nt) = P1dt(:);
    udt(4*Nt+1 : 5*Nt) = P2dt(:);


end
