function udt = eom2nr(t, u, params)
    % non-relativistic, small-gradient approximation to the equations of motion
    
%   choice of Pi: 
%   Pi_xx =  2 P1 (1+UX^2) +  2 UX UY P2
%   Pi_yy = -2 P1 (1+UY^2) +  2 UX UY P2
%   Pi_xy =  P2 ( 2 + UX^2 + UY^2 )

    A = params.A;
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

    m11 = (-3).*T.^2;
    m12 = (-2).*P1.*ux+(-3/2).*T.^3.*ux+(-2).*P2.*uy;
    m13 = (-2).*P2.*ux+2.*P1.*uy+(-3/2).*T.^3.*uy;


    m22 = 2.*P1+(3/2).*T.^3;
    m23 = 2.*P2;
    m24 = 2.*ux;
    m25 = 2.*uy;


    m32 = 2.*P2;
    m33 = (-2).*P1+(3/2).*T.^3;
    m34 = (-2).*uy;
    m35 = 2.*ux;

    m42 = (1/6).*T.^(-4).*(16.*(1+A).*P1.^2.*ux+3.*T.^3.*(T.^3.*ux+4.*P2.*uy)+2.*(1+A).*P1.*(9.*T.^3.*ux+8.*P2.*uy));
    m43 = (1/6).*T.^(-4).*((-16).*(1+A).*P1.^2.*uy+2.*(1+A).*P1.*(8.*P2.*ux+9.*T.^3.*uy)+(-3).*(4.*P2.*T.^3.*ux+T.^6.*uy));
    m44 = 2.*(1+A).*T.^(-1);

    m52 = (1/6).*T.^(-4).*(18.*(1+A).*P2.*T.^3.*ux+16.*(1+A).*P2.^2.*uy+3.*T.^6.*uy+4.*P1.*(4.*(1+A).*P2.*ux+(-3).*T.^3.*uy));
    m53 = (1/6).*T.^(-4).*(16.*(1+A).*P2.^2.*ux+3.*T.^3.*(4.*P1+T.^3).*ux+(-2).*(1+A).*P2.*(8.*P1+(-9).*T.^3).*uy);
    m55 = 2.*(1+A).*T.^(-1);


    R1 = (1/2).*((-6).*T.^2.*(Tdx.*ux+Tdy.*uy)+(-4).*(P2.*(uxdy+uydx)+P1.*(uxdx+(-1).*uydy))+(-3).*T.^3.*(uxdx+uydy));
    R2 = (1/4).*(8.*P1dx+8.*P2dy+6.*T.^2.*Tdx+4.*P1.*ux.*uxdx+3.*T.^3.*ux.*uxdx+(-4).*P2.*ux.*uxdy+8.*P2.*uxdx.*uy+6.*T.^3.*uxdy.*uy+(-4).*P2.*ux.*uydx+12.*P1.*ux.*uydy+(-3).*T.^3.*ux.*uydy+8.*P2.*uy.*uydy);
    R3 = (1/4).*((-8).*P1dy+8.*P2dx+6.*T.^2.*Tdy+8.*P2.*ux.*uxdx+(-12).*P1.*uxdx.*uy+(-3).*T.^3.*uxdx.*uy+(-4).*P2.*uxdy.*uy+6.*T.^3.*ux.*uydx+(-4).*P2.*uy.*uydx+8.*P2.*ux.*uydy+(-4).*P1.*uy.*uydy+3.*T.^3.*uy.*uydy);
    R4 = (1/6).*T.^(-4).*(3.*T.^3.*(4.*(1+A).*P1dx.*ux+4.*((-1).*A.*P2.*uxdy+P1dy.*uy+A.*P1dy.*uy+A.*P2.*uydx)+T.^3.*(uxdx+(-1).*uydy))+16.*(1+A).*P1.^2.*(uxdx+(-1).*uydy)+2.*P1.*(6.*T.^4+8.*(1+A).*P2.*(uxdy+uydx)+9.*(1+A).*T.^3.*(uxdx+uydy)));
    R5 = (1/6).*T.^(-4).*(16.*(1+A).*P2.^2.*(uxdy+uydx)+3.*T.^3.*(4.*(1+A).*P2dx.*ux+T.^3.*uxdy+4.*P2dy.*uy+T.^3.*uydx+4.*A.*(P1.*uxdy+P2dy.*uy+(-1).*P1.*uydx))+2.*P2.*(6.*T.^4+8.*(1+A).*P1.*(uxdx+(-1).*uydy)+9.*(1+A).*T.^3.*(uxdx+uydy)));

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

    i11 = i11./mdet;
    i12 = i12./mdet;
    i13 = i13./mdet;
    i14 = i14./mdet;
    i15 = i15./mdet;

    i22 = i22./mdet;
    i23 = i23./mdet;
    i24 = i24./mdet;
    i25 = i25./mdet;

    i32 = i32./mdet;
    i33 = i33./mdet;
    i34 = i34./mdet;
    i35 = i35./mdet;

    i42 = i42./mdet;
    i43 = i43./mdet;
    i44 = i44./mdet;
    i45 = i45./mdet;

    i52 = i52./mdet;
    i53 = i53./mdet;
    i54 = i54./mdet;
    i55 = i55./mdet;

    
    Tdt  = -(i11.*R1+i12.*R2+i13.*R3+i14.*R4+i15.*R5);
    uxdt = -(i22.*R2+i23.*R3+i24.*R4+i25.*R5);
    uydt = -(i32.*R2+i33.*R3+i34.*R4+i35.*R5);
    P1dt = -(i42.*R2+i43.*R3+i44.*R4+i45.*R5);
    P2dt = -(i52.*R2+i53.*R3+i54.*R4+i55.*R5);


    udt = zeros(Nx.*Ny.*5, 1);

    udt(0*Nt+1 : 1*Nt) = Tdt(:);
    udt(1*Nt+1 : 2*Nt) = uxdt(:);
    udt(2*Nt+1 : 3*Nt) = uydt(:);
    udt(3*Nt+1 : 4*Nt) = P1dt(:);
    udt(4*Nt+1 : 5*Nt) = P2dt(:);

    % udt(0*Nt+1 : 1*Nt) = gather(Tdt(:));
    % udt(1*Nt+1 : 2*Nt) = gather(uxdt(:));
    % udt(2*Nt+1 : 3*Nt) = gather(uydt(:));
    % udt(3*Nt+1 : 4*Nt) = gather(P1dt(:));
    % udt(4*Nt+1 : 5*Nt) = gather(P2dt(:));


end
