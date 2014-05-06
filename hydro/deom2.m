function udtdt = deom2(t, u, udt, params)
%   choice of Pi: 
%   Pi_xx =  2 P1 (1+UX^2) +  2 UX UY P2
%   Pi_yy = -2 P1 (1+UY^2) +  2 UX UY P2
%   Pi_xy =  P2 ( 2 + UX^2 + UY^2 )

    A = params.A;
    p1 = A+1;
    p2 = 2.*A;
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


    % time derivatives
    Tdt  = reshape(udt(0*Nt+1 : 1*Nt), Nx, Ny);
    uxdt = reshape(udt(1*Nt+1 : 2*Nt), Nx, Ny);
    uydt = reshape(udt(2*Nt+1 : 3*Nt), Nx, Ny);
    P1dt = reshape(udt(3*Nt+1 : 4*Nt), Nx, Ny);
    P2dt = reshape(udt(4*Nt+1 : 5*Nt), Nx, Ny);

    utdt = (uxdt.*ux + uydt.*uy)./ut;

    % compute space dervatives
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

    utdx = (1/2).*ut.^(-1).*(2.*ux.*uxdx+2.*uy.*uydx);
    utdy = (1/2).*ut.^(-1).*(2.*ux.*uxdy+2.*uy.*uydy);

    % compute space derivatives of dt
    Tdtdx  = dx(Tdt);
    Tdtdy  = dy(Tdt);
    uxdtdx = dx(uxdt);
    uxdtdy = dy(uxdt);
    uydtdx = dx(uydt);
    uydtdy = dy(uydt);
    P1dtdx = dx(P1dt);
    P1dtdy = dy(P1dt);
    P2dtdx = dx(P2dt);
    P2dtdy = dy(P2dt);

    utdtdx = dx(utdt);
    utdtdy = dy(utdt);



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



    global xt

    xt.tic('m');

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




    mdt11 = (-6).*T.*Tdt.*ut2+(-6).*T2.*ut.*utdt;
    mdt12 = (-9/2).*T2.*Tdt.*ux+(-3/2).*T3.*uxdt+P2dt.*uy.*((-2)+ux2+(-1).*uy2)+(-2).*P1dt.*ux.*(1+uy2)+(-2).*P1.*uxdt.*(1+uy2)+(-4).*P1.*ux.*uy.*uydt+P2.*((-2)+ux2+(-1).*uy2).*uydt+P2.*uy.*(2.*ux.*uxdt+(-2).*uy.*uydt);
    mdt13 = (1/2).*((-9).*T2.*Tdt+4.*P1dt.*(1+ux2)+8.*P1.*ux.*uxdt).*uy+P2dt.*ux.*((-2)+(-1).*ux2+uy2)+P2.*uxdt.*((-2)+(-1).*ux2+uy2)+(1/2).*((-3).*T3+4.*P1.*(1+ux2)).*uydt+P2.*ux.*((-2).*ux.*uxdt+2.*uy.*uydt);

    mdt22 = (-1/2).*ut.^(-3).*utdt.*ux.*uy.*(3.*T3.*ut2.*(2+ux2+2.*uy2)+2.*P2.*ux.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+4.*P1.*(ux4.*(1+(-3).*uy2)+2.*(1+uy2)+(-3).*ux2.*((-1)+uy4)))+(1/4).*ut.^(-2).*uxdt.*uy.*(3.*T3.*ut2.*(2+ux2+2.*uy2)+2.*P2.*ux.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+4.*P1.*(ux4.*(1+(-3).*uy2)+2.*(1+uy2)+(-3).*ux2.*((-1)+uy4)))+(1/4).*ut.^(-2).*ux.*(3.*T3.*ut2.*(2+ux2+2.*uy2)+2.*P2.*ux.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+4.*P1.*(ux4.*(1+(-3).*uy2)+2.*(1+uy2)+(-3).*ux2.*((-1)+uy4))).*uydt+(1/4).*ut.^(-2).*ux.*uy.*(9.*T2.*Tdt.*ut2.*(2+ux2+2.*uy2)+6.*T3.*ut.*utdt.*(2+ux2+2.*uy2)+2.*P2dt.*ux.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+2.*P2.*uxdt.*uy.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4)+4.*P1dt.*(ux4.*(1+(-3).*uy2)+2.*(1+uy2)+(-3).*ux2.*((-1)+uy4))+2.*P2.*ux.*(2+3.*(ux2+ux4)+uy2+(-3).*uy4).*uydt+3.*T3.*ut2.*(2.*ux.*uxdt+4.*uy.*uydt)+2.*P2.*ux.*uy.*(3.*(2.*ux.*uxdt+4.*ux3.*uxdt)+2.*uy.*uydt+(-12).*uy3.*uydt)+4.*P1.*(4.*ux3.*uxdt.*(1+(-3).*uy2)+(-6).*ux.*uxdt.*((-1)+uy4)+4.*uy.*uydt+(-6).*ux4.*uy.*uydt+(-12).*ux2.*uy3.*uydt));
    mdt23 = (-1/2).*ut.^(-3).*utdt.*ux.*uy.*(ux.*uy.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2))+2.*P2.*(4+(-3).*ux4+(-3).*ux.^6+6.*uy2+4.*uy4+ux2.*(4+3.*(uy2+uy4))))+(1/4).*ut.^(-2).*uxdt.*uy.*(ux.*uy.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2))+2.*P2.*(4+(-3).*ux4+(-3).*ux.^6+6.*uy2+4.*uy4+ux2.*(4+3.*(uy2+uy4))))+(1/4).*ut.^(-2).*ux.*(ux.*uy.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2))+2.*P2.*(4+(-3).*ux4+(-3).*ux.^6+6.*uy2+4.*uy4+ux2.*(4+3.*(uy2+uy4)))).*uydt+(1/4).*ut.^(-2).*ux.*uy.*(uxdt.*uy.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2))+2.*P2dt.*(4+(-3).*ux4+(-3).*ux.^6+6.*uy2+4.*uy4+ux2.*(4+3.*(uy2+uy4)))+ux.*((-3).*T3.*ut2+4.*P1.*(1+ux2).*(1+3.*ux2+3.*uy2)).*uydt+ux.*uy.*((-9).*T2.*Tdt.*ut2+(-6).*T3.*ut.*utdt+4.*P1dt.*(1+ux2).*(1+3.*ux2+3.*uy2)+8.*P1.*ux.*uxdt.*(1+3.*ux2+3.*uy2)+4.*P1.*(1+ux2).*(6.*ux.*uxdt+6.*uy.*uydt))+2.*P2.*((-12).*ux3.*uxdt+(-18).*ux.^5.*uxdt+2.*ux.*uxdt.*(4+3.*(uy2+uy4))+12.*uy.*uydt+16.*uy3.*uydt+3.*ux2.*(2.*uy.*uydt+4.*uy3.*uydt)));
    mdt24 = 4.*ux3.*uxdt.*uy+4.*ux.*(1+ux2).*uxdt.*uy+2.*ux2.*(1+ux2).*uydt;
    mdt25 = uxdt.*uy2.*(2+3.*ux2+uy2)+2.*ux.*uy.*(2+3.*ux2+uy2).*uydt+ux.*uy2.*(6.*ux.*uxdt+2.*uy.*uydt);

    mdt32 = (-1/2).*ut.^(-3).*utdt.*ux.*uy.*((-1).*ux.*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2))+2.*P2.*(3.*ux2.*(2+uy2)+ux4.*(4+3.*uy2)+(-1).*(1+uy2).*((-4)+3.*uy4)))+(1/4).*ut.^(-2).*uxdt.*uy.*((-1).*ux.*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2))+2.*P2.*(3.*ux2.*(2+uy2)+ux4.*(4+3.*uy2)+(-1).*(1+uy2).*((-4)+3.*uy4)))+(1/4).*ut.^(-2).*ux.*((-1).*ux.*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2))+2.*P2.*(3.*ux2.*(2+uy2)+ux4.*(4+3.*uy2)+(-1).*(1+uy2).*((-4)+3.*uy4))).*uydt+(1/4).*ut.^(-2).*ux.*uy.*((-1).*uxdt.*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2))+2.*P2dt.*(3.*ux2.*(2+uy2)+ux4.*(4+3.*uy2)+(-1).*(1+uy2).*((-4)+3.*uy4))+(-1).*ux.*(3.*T3.*ut2+4.*P1.*(1+uy2).*(1+3.*ux2+3.*uy2)).*uydt+2.*P2.*(6.*ux.*uxdt.*(2+uy2)+4.*ux3.*uxdt.*(4+3.*uy2)+6.*ux2.*uy.*uydt+6.*ux4.*uy.*uydt+(-12).*uy3.*(1+uy2).*uydt+(-2).*uy.*((-4)+3.*uy4).*uydt)+(-1).*ux.*uy.*(9.*T2.*Tdt.*ut2+6.*T3.*ut.*utdt+4.*P1dt.*(1+uy2).*(1+3.*ux2+3.*uy2)+8.*P1.*uy.*(1+3.*ux2+3.*uy2).*uydt+4.*P1.*(1+uy2).*(6.*ux.*uxdt+6.*uy.*uydt)));
    mdt33 = (-1/2).*ut.^(-3).*utdt.*ux.*uy.*(3.*T3.*ut2.*(2+2.*ux2+uy2)+4.*P1.*((-2).*(1+ux2)+3.*((-1)+ux4).*uy2+((-1)+3.*ux2).*uy4)+2.*P2.*ux.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4)))+(1/4).*ut.^(-2).*uxdt.*uy.*(3.*T3.*ut2.*(2+2.*ux2+uy2)+4.*P1.*((-2).*(1+ux2)+3.*((-1)+ux4).*uy2+((-1)+3.*ux2).*uy4)+2.*P2.*ux.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4)))+(1/4).*ut.^(-2).*ux.*(3.*T3.*ut2.*(2+2.*ux2+uy2)+4.*P1.*((-2).*(1+ux2)+3.*((-1)+ux4).*uy2+((-1)+3.*ux2).*uy4)+2.*P2.*ux.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4))).*uydt+(1/4).*ut.^(-2).*ux.*uy.*(9.*T2.*Tdt.*ut2.*(2+2.*ux2+uy2)+6.*T3.*ut.*utdt.*(2+2.*ux2+uy2)+4.*P1dt.*((-2).*(1+ux2)+3.*((-1)+ux4).*uy2+((-1)+3.*ux2).*uy4)+2.*P2dt.*ux.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4))+2.*P2.*uxdt.*uy.*(2+ux2+(-3).*ux4+3.*(uy2+uy4))+2.*P2.*ux.*(2+ux2+(-3).*ux4+3.*(uy2+uy4)).*uydt+3.*T3.*ut2.*(4.*ux.*uxdt+2.*uy.*uydt)+4.*P1.*((-4).*ux.*uxdt+12.*ux3.*uxdt.*uy2+6.*ux.*uxdt.*uy4+6.*((-1)+ux4).*uy.*uydt+4.*((-1)+3.*ux2).*uy3.*uydt)+2.*P2.*ux.*uy.*(2.*ux.*uxdt+(-12).*ux3.*uxdt+3.*(2.*uy.*uydt+4.*uy3.*uydt)));
    mdt34 = (-2).*uxdt.*uy2.*(1+uy2)+(-4).*ux.*uy3.*uydt+(-4).*ux.*uy.*(1+uy2).*uydt;
    mdt35 = 2.*ux.*uxdt.*uy.*(2+ux2+3.*uy2)+ux2.*(2+ux2+3.*uy2).*uydt+ux2.*uy.*(2.*ux.*uxdt+6.*uy.*uydt);

    mdt42 = (1/12).*(9.*T2.*Tdt.*(T3.*ux.*(2+ux2+3.*uy2)+A2.*uy.*(2.*P1.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-1).*P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)))+3.*T3.*(3.*T2.*Tdt.*ux.*(2+ux2+3.*uy2)+T3.*uxdt.*(2+ux2+3.*uy2)+A2.*(2.*P1.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-1).*P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)).*uydt+T3.*ux.*(2.*ux.*uxdt+6.*uy.*uydt)+A2.*uy.*(2.*P1dt.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+2.*P1.*uxdt.*(ux+(-1).*uy).*uy.*(ux+uy)+(-1).*P2dt.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)+2.*P1.*ux.*uy.*(ux+uy).*(uxdt+(-1).*uydt)+2.*P1.*ux.*(ux+(-1).*uy).*(ux+uy).*uydt+2.*P1.*ux.*(ux+(-1).*uy).*uy.*(uxdt+uydt)+(-1).*P2.*(4.*ux.*(2+ux2).*uxdt+(-4).*ux.*uxdt.*uy2+(-4).*((-2)+ux2).*uy.*uydt+4.*uy3.*uydt)))+2.*A1.*(16.*P1.*P1dt.*ux.*(1+uy2).*(2+ux2+uy2)+8.*P1.^2.*uxdt.*(1+uy2).*(2+ux2+uy2)+6.*P2dt.*T3.*uy.*(2+ux4+3.*uy2+uy4+ux2.*(1+(-2).*uy2))+18.*P2.*T2.*Tdt.*uy.*(2+ux4+3.*uy2+uy4+ux2.*(1+(-2).*uy2))+P1dt.*(3.*T3.*ux.*(6+7.*uy2+4.*uy4+ux2.*(3+(-4).*uy2))+4.*P2.*uy.*((-1).*ux4+(2+uy2).^2))+16.*P1.^2.*ux.*uy.*(2+ux2+uy2).*uydt+6.*P2.*T3.*(2+ux4+3.*uy2+uy4+ux2.*(1+(-2).*uy2)).*uydt+8.*P1.^2.*ux.*(1+uy2).*(2.*ux.*uxdt+2.*uy.*uydt)+6.*P2.*T3.*uy.*(4.*ux3.*uxdt+2.*ux.*uxdt.*(1+(-2).*uy2)+6.*uy.*uydt+(-4).*ux2.*uy.*uydt+4.*uy3.*uydt)+P1.*(9.*T2.*Tdt.*ux.*(6+7.*uy2+4.*uy4+ux2.*(3+(-4).*uy2))+3.*T3.*uxdt.*(6+7.*uy2+4.*uy4+ux2.*(3+(-4).*uy2))+4.*P2dt.*uy.*((-1).*ux4+(2+uy2).^2)+4.*P2.*((-1).*ux4+(2+uy2).^2).*uydt+3.*T3.*ux.*(2.*ux.*uxdt.*(3+(-4).*uy2)+14.*uy.*uydt+(-8).*ux2.*uy.*uydt+16.*uy3.*uydt)+4.*P2.*uy.*((-4).*ux3.*uxdt+4.*uy.*(2+uy2).*uydt))));
    mdt43 = (1/12).*(9.*T2.*Tdt.*((-1).*T3.*uy.*(2+3.*ux2+uy2)+A2.*ux.*(2.*P1.*ux.*uy.*((-1).*ux2+uy2)+P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)))+3.*T3.*((-3).*T2.*Tdt.*uy.*(2+3.*ux2+uy2)+A2.*uxdt.*(2.*P1.*ux.*uy.*((-1).*ux2+uy2)+P2.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4))+(-1).*T3.*(2+3.*ux2+uy2).*uydt+(-1).*T3.*uy.*(6.*ux.*uxdt+2.*uy.*uydt)+A2.*ux.*(2.*P1dt.*ux.*uy.*((-1).*ux2+uy2)+2.*P1.*uxdt.*uy.*((-1).*ux2+uy2)+P2dt.*((2+ux2).^2+(-2).*((-2)+ux2).*uy2+uy4)+2.*P1.*ux.*((-1).*ux2+uy2).*uydt+2.*P1.*ux.*uy.*((-2).*ux.*uxdt+2.*uy.*uydt)+P2.*(4.*ux.*(2+ux2).*uxdt+(-4).*ux.*uxdt.*uy2+(-4).*((-2)+ux2).*uy.*uydt+4.*uy3.*uydt)))+2.*A1.*((-16).*P1.*P1dt.*(1+ux2).*uy.*(2+ux2+uy2)+(-16).*P1.^2.*ux.*uxdt.*uy.*(2+ux2+uy2)+(-6).*P2dt.*T3.*ux.*(2+ux4+uy2+uy4+ux2.*(3+(-2).*uy2))+(-18).*P2.*T2.*Tdt.*ux.*(2+ux4+uy2+uy4+ux2.*(3+(-2).*uy2))+(-6).*P2.*T3.*uxdt.*(2+ux4+uy2+uy4+ux2.*(3+(-2).*uy2))+P1dt.*(4.*P2.*ux.*((2+ux2).^2+(-1).*uy4)+3.*T3.*uy.*(ux2.*(7+4.*ux2+(-4).*uy2)+3.*(2+uy2)))+(-8).*P1.^2.*(1+ux2).*(2+ux2+uy2).*uydt+(-8).*P1.^2.*(1+ux2).*uy.*(2.*ux.*uxdt+2.*uy.*uydt)+(-6).*P2.*T3.*ux.*(4.*ux3.*uxdt+2.*ux.*uxdt.*(3+(-2).*uy2)+2.*uy.*uydt+(-4).*ux2.*uy.*uydt+4.*uy3.*uydt)+P1.*(4.*P2dt.*ux.*((2+ux2).^2+(-1).*uy4)+4.*P2.*uxdt.*((2+ux2).^2+(-1).*uy4)+9.*T2.*Tdt.*uy.*(ux2.*(7+4.*ux2+(-4).*uy2)+3.*(2+uy2))+3.*T3.*(ux2.*(7+4.*ux2+(-4).*uy2)+3.*(2+uy2)).*uydt+4.*P2.*ux.*(4.*ux.*(2+ux2).*uxdt+(-4).*uy3.*uydt)+3.*T3.*uy.*(2.*ux.*uxdt.*(7+4.*ux2+(-4).*uy2)+6.*uy.*uydt+ux2.*(8.*ux.*uxdt+(-8).*uy.*uydt)))));
    mdt44 = 3.*A1.*T2.*Tdt.*ut2.*(2+ux2+uy2)+2.*A1.*T3.*ut.*utdt.*(2+ux2+uy2)+A1.*T3.*ut2.*(2.*ux.*uxdt+2.*uy.*uydt);

    mdt52 = (1/6).*uxdt.*uy.*(3.*T3.*uy.*(T3.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)))+A1.*(4.*P1.*(1+uy2).*((-3).*T3.*(1+2.*ux2).*uy+2.*P2.*ux.*(2+ux2+uy2))+P2.*(4.*P2.*uy.*((-1).*ux4+(2+uy2).^2)+3.*T3.*ux.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2)))))+(1/6).*ux.*(3.*T3.*uy.*(T3.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)))+A1.*(4.*P1.*(1+uy2).*((-3).*T3.*(1+2.*ux2).*uy+2.*P2.*ux.*(2+ux2+uy2))+P2.*(4.*P2.*uy.*((-1).*ux4+(2+uy2).^2)+3.*T3.*ux.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2))))).*uydt+(1/6).*ux.*uy.*(9.*T2.*Tdt.*uy.*(T3.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)))+3.*T3.*(T3.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2))).*uydt+3.*T3.*uy.*(3.*T2.*Tdt.*(1+uy2)+2.*T3.*uy.*uydt+A2.*(2.*P1dt.*(1+ux2).*(1+uy2)+4.*P1.*ux.*uxdt.*(1+uy2)+P2dt.*ux.*uy.*((-1).*ux2+uy2)+P2.*uxdt.*uy.*((-1).*ux2+uy2)+4.*P1.*(1+ux2).*uy.*uydt+P2.*ux.*((-1).*ux2+uy2).*uydt+P2.*ux.*uy.*((-2).*ux.*uxdt+2.*uy.*uydt)))+A1.*(4.*P1dt.*(1+uy2).*((-3).*T3.*(1+2.*ux2).*uy+2.*P2.*ux.*(2+ux2+uy2))+P2dt.*(4.*P2.*uy.*((-1).*ux4+(2+uy2).^2)+3.*T3.*ux.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2)))+8.*P1.*uy.*((-3).*T3.*(1+2.*ux2).*uy+2.*P2.*ux.*(2+ux2+uy2)).*uydt+4.*P1.*(1+uy2).*((-9).*T2.*Tdt.*(1+2.*ux2).*uy+(-12).*T3.*ux.*uxdt.*uy+2.*P2dt.*ux.*(2+ux2+uy2)+2.*P2.*uxdt.*(2+ux2+uy2)+(-3).*T3.*(1+2.*ux2).*uydt+2.*P2.*ux.*(2.*ux.*uxdt+2.*uy.*uydt))+P2.*(4.*P2dt.*uy.*((-1).*ux4+(2+uy2).^2)+9.*T2.*Tdt.*ux.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2))+3.*T3.*uxdt.*(6+uy2+(-4).*uy4+ux2.*(5+4.*uy2))+4.*P2.*((-1).*ux4+(2+uy2).^2).*uydt+3.*T3.*ux.*(2.*ux.*uxdt.*(5+4.*uy2)+2.*uy.*uydt+8.*ux2.*uy.*uydt+(-16).*uy3.*uydt)+4.*P2.*uy.*((-4).*ux3.*uxdt+4.*uy.*(2+uy2).*uydt))));
    mdt53 = (1/6).*uxdt.*uy.*(3.*T3.*ux.*(T3.*(1+ux2)+A2.*(P2.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1.*(1+ux2).*(1+uy2)))+A1.*(12.*P1.*T3.*ux.*(1+ux2).*(1+2.*uy2)+4.*P2.^2.*ux.*((2+ux2).^2+(-1).*uy4)+P2.*uy.*((-8).*P1.*(1+ux2).*(2+ux2+uy2)+3.*T3.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2))))+(1/6).*ux.*(3.*T3.*ux.*(T3.*(1+ux2)+A2.*(P2.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1.*(1+ux2).*(1+uy2)))+A1.*(12.*P1.*T3.*ux.*(1+ux2).*(1+2.*uy2)+4.*P2.^2.*ux.*((2+ux2).^2+(-1).*uy4)+P2.*uy.*((-8).*P1.*(1+ux2).*(2+ux2+uy2)+3.*T3.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2)))).*uydt+(1/6).*ux.*uy.*(9.*T2.*Tdt.*ux.*(T3.*(1+ux2)+A2.*(P2.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1.*(1+ux2).*(1+uy2)))+3.*T3.*uxdt.*(T3.*(1+ux2)+A2.*(P2.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1.*(1+ux2).*(1+uy2)))+3.*T3.*ux.*(3.*T2.*Tdt.*(1+ux2)+2.*T3.*ux.*uxdt+A2.*(P2dt.*ux.*(ux+(-1).*uy).*uy.*(ux+uy)+P2.*uxdt.*(ux+(-1).*uy).*uy.*(ux+uy)+(-2).*P1dt.*(1+ux2).*(1+uy2)+(-4).*P1.*ux.*uxdt.*(1+uy2)+P2.*ux.*uy.*(ux+uy).*(uxdt+(-1).*uydt)+(-4).*P1.*(1+ux2).*uy.*uydt+P2.*ux.*(ux+(-1).*uy).*(ux+uy).*uydt+P2.*ux.*(ux+(-1).*uy).*uy.*(uxdt+uydt)))+A1.*(12.*P1dt.*T3.*ux.*(1+ux2).*(1+2.*uy2)+36.*P1.*T2.*Tdt.*ux.*(1+ux2).*(1+2.*uy2)+24.*P1.*T3.*ux2.*uxdt.*(1+2.*uy2)+12.*P1.*T3.*(1+ux2).*uxdt.*(1+2.*uy2)+8.*P2.*P2dt.*ux.*((2+ux2).^2+(-1).*uy4)+4.*P2.^2.*uxdt.*((2+ux2).^2+(-1).*uy4)+P2dt.*uy.*((-8).*P1.*(1+ux2).*(2+ux2+uy2)+3.*T3.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2))+48.*P1.*T3.*ux.*(1+ux2).*uy.*uydt+P2.*((-8).*P1.*(1+ux2).*(2+ux2+uy2)+3.*T3.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2)).*uydt+4.*P2.^2.*ux.*(4.*ux.*(2+ux2).*uxdt+(-4).*uy3.*uydt)+P2.*uy.*((-8).*P1dt.*(1+ux2).*(2+ux2+uy2)+(-16).*P1.*ux.*uxdt.*(2+ux2+uy2)+9.*T2.*Tdt.*(6+ux2+(-4).*ux4+(5+4.*ux2).*uy2)+(-8).*P1.*(1+ux2).*(2.*ux.*uxdt+2.*uy.*uydt)+3.*T3.*(2.*ux.*uxdt+(-16).*ux3.*uxdt+8.*ux.*uxdt.*uy2+2.*(5+4.*ux2).*uy.*uydt))));
    mdt55 = 3.*A1.*T2.*Tdt.*ut2.*ux.*uy.*(2+ux2+uy2)+2.*A1.*T3.*ut.*utdt.*ux.*uy.*(2+ux2+uy2)+A1.*T3.*ut2.*uxdt.*uy.*(2+ux2+uy2)+A1.*T3.*ut2.*ux.*(2+ux2+uy2).*uydt+A1.*T3.*ut2.*ux.*uy.*(2.*ux.*uxdt+2.*uy.*uydt);


    xt.toc('m');
 




    xt.tic('rp');

    RdtTdt = (-3/2).*T.*ut.*(4.*Tdx.*ux+4.*Tdy.*uy+3.*T.*(uxdx+ut.*(utdx.*ux+utdy.*uy)+(-1).*ux.*(ux.*uxdx+uy.*(uxdy+uydx))+uydy+(-1).*uy2.*uydy));
    RdtTdtdx = (-3).*T2.*ut.*ux;
    RdtTdtdy = (-3).*T2.*ut.*uy;
    Rdtutdt = (1/2).*((-6).*T2.*(Tdx.*ux+Tdy.*uy)+(-3).*T3.*(uxdx+2.*ut.*(utdx.*ux+utdy.*uy)+(-1).*ux.*(ux.*uxdx+uy.*(uxdy+uydx))+uydy+(-1).*uy2.*uydy)+(-4).*P1.*((1+ux2).*uxdx+(-1).*(1+uy2).*uydy)+(-2).*P2.*(4.*uxdy+4.*uydx+3.*ux2.*(uxdy+uydx)+3.*uy2.*(uxdy+uydx)+2.*ux.*uy.*(uxdx+uydy)));
    Rdtutdtdx = (-3/2).*T3.*ut2.*ux+2.*P1.*(ux+ux3)+P2.*uy.*(2+3.*ux2+uy2);
    Rdtutdtdy = P2.*ux.*(2+ux2+3.*uy2)+(-1/2).*uy.*(3.*T3.*ut2+4.*P1.*(1+uy2));
    Rdtuxdt = (1/2).*((-6).*T2.*Tdx.*ut+4.*P1.*(utdx+3.*utdx.*ux2+(-2).*ut.*ux.*uxdx)+3.*T3.*ut.*((-1).*ut.*utdx+2.*ux.*uxdx+uy.*(uxdy+uydx))+2.*P2.*(6.*utdx.*ux.*uy+utdy.*(2+3.*ux2+3.*uy2)+(-2).*ut.*uy.*(uxdx+uydy)));
    Rdtuxdtdx = (-1/2).*ut.*((-3).*T3.*((-1)+ux2)+4.*P1.*(1+ux2)+4.*P2.*ux.*uy);
    Rdtuxdtdy = (1/2).*ut.*((-2).*P2.*(1+ut2)+3.*T3.*ux.*uy);
    Rdtuydt = (1/2).*((-6).*T2.*Tdy.*ut+3.*T3.*ut.*((-1).*ut.*utdy+ux.*(uxdy+uydx)+2.*uy.*uydy)+(-4).*P1.*(utdy+3.*utdy.*uy2+(-2).*ut.*uy.*uydy)+2.*P2.*(6.*utdy.*ux.*uy+utdx.*(2+3.*ux2+3.*uy2)+(-2).*ut.*ux.*(uxdx+uydy)));
    Rdtuydtdx = (1/2).*ut.*((-2).*P2.*(1+ut2)+3.*T3.*ux.*uy);
    Rdtuydtdy = (1/2).*ut.*((-4).*P2.*ux.*uy+3.*T3.*((-1)+uy2)+4.*P1.*(1+uy2));
    RdtP1dt = 2.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy));
    RdtP1dtdx = 0;
    RdtP1dtdy = 0;
    RdtP2dt = utdx.*uy.*(2+3.*ux2+uy2)+utdy.*ux.*(2+ux2+3.*uy2)+(-1).*ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+2.*ux.*uy.*(uxdx+uydy));
    RdtP2dtdx = 0;
    RdtP2dtdy = 0;

    Rdt1 = RdtTdt.*Tdt + RdtTdtdx.*Tdtdx + RdtTdtdy.*Tdtdy + Rdtutdt.*utdt + Rdtutdtdx.*utdtdx + Rdtutdtdy.*utdtdy + Rdtuxdt.*uxdt + Rdtuxdtdx.*uxdtdx + Rdtuxdtdy.*uxdtdy + Rdtuydt.*uydt + Rdtuydtdx.*uydtdx + Rdtuydtdy.*uydtdy + RdtP1dt.*P1dt + RdtP1dtdx.*P1dtdx + RdtP1dtdy.*P1dtdy + RdtP2dt.*P2dt + RdtP2dtdx.*P2dtdx + RdtP2dtdy.*P2dtdy;


    RdtTdt = (3/4).*T.*ux.*uy.*(2.*ut.*(2.*Tdx+3.*T.*uxdy.*uy)+(-3).*T.*ux.*(3.*ut2.*utdx.*ux+3.*ut2.*utdy.*uy+ut.*((-1).*uxdx+(-3).*ux.*(ux.*uxdx+uy.*(uxdy+uydx))+uydy+(-3).*uy2.*uydy)));
    RdtTdtdx = (3/2).*T2.*ut.*ux.*uy;
    RdtTdtdy = 0;
    Rdtutdt = (1/4).*(8.*P1.*uy3.*uydy+ux.*(uy.*(8.*P1dx+6.*T2.*Tdx+P2dy.*(16+12.*uy2)+uy.*(8.*P1.*uxdy+6.*T3.*uxdy+P2.*(8.*uxdx+(-8).*uydy)))+ux.*(uy.*(2.*(4.*P2dx.*uy+P2.*((-10).*uxdy+(-2).*uydx+(-3).*uy2.*(uxdy+uydx)))+4.*P1.*(uxdx+uydy+3.*uy2.*uydy)+3.*T3.*(uxdx+(-6).*ut.*utdy.*uy+((-1)+3.*uy2).*uydy))+ux.*((-8).*P1.*uxdy+3.*ux.*uy.*((-4).*P1.*uxdx+3.*T3.*uxdx+(-2).*P2.*(uxdy+uydx))+uy.*(8.*P1dx+12.*P2dy+(-18).*T3.*ut.*utdx+3.*uy.*(3.*T3.*(uxdy+uydx)+(-4).*P2.*(uxdx+uydy)))))));
    Rdtutdtdx = (3/4).*ux2.*uy.*((-3).*T3.*ut2.*ux+4.*P1.*(ux+ux3)+2.*P2.*uy.*(2+3.*ux2+uy2));
    Rdtutdtdy = (1/4).*(P1.*(8.*ux4+(-8).*uy2.*(1+uy2)+ux2.*(8+(-12).*uy2+(-12).*uy4))+ux.*uy.*((-9).*T3.*ut2.*ux.*uy+P2.*(16+16.*uy2+ux2.*(28+6.*ux2+18.*uy2))));
    Rdtuxdt = (1/4).*(2.*uy.*(4.*(P1dx+P2dy).*ut+3.*T2.*Tdx.*ut+8.*P2.*utdy+uy.*(3.*T3.*ut.*uxdy+2.*(P2dy.*ut+4.*P2.*utdy).*uy+4.*ut.*(P1.*uxdy+P2.*(uxdx+(-1).*uydy))))+ux.*(16.*P1.*utdy+ux.*(P1.*((-24).*ut.*uxdy+36.*utdx.*uy+ux.*(32.*utdy+60.*utdx.*ux.*uy+(-48).*ut.*uxdx.*uy))+uy.*(4.*(6.*P1dx+P2dy).*ut+(-27).*T3.*utdx+84.*P2.*utdy+ux.*(T3.*((-27).*utdx.*ux+36.*ut.*uxdx+(-18).*utdy.*uy)+6.*P2.*(5.*utdy.*ux+12.*utdx.*uy+(-4).*ut.*(uxdy+uydx)))+9.*uy.*((-3).*(T3.*utdx+(-2).*P2.*utdy).*uy+3.*T3.*ut.*(uxdy+uydx)+(-4).*P2.*ut.*(uxdx+uydy))))+2.*uy.*(3.*T3.*ut.*(uxdx+(-1).*uydy)+4.*ut.*((-1).*P2.*(5.*uxdy+uydx)+P1.*(uxdx+uydy))+uy.*(8.*P2dx.*ut+3.*(2.*P2.*(utdx.*(2+uy2)+(-1).*ut.*uy.*(uxdy+uydx))+(-1).*(4.*P1+3.*T3).*(utdy+utdy.*uy2+(-1).*ut.*uy.*uydy))))));
    Rdtuxdtdx = (-1/4).*ut.*ux.*uy.*(P1.*ux.*((-4)+12.*ux2)+(-3).*T3.*(ux+3.*ux3)+4.*P2.*((-2)+3.*ux2).*uy);
    Rdtuxdtdy = (-1/4).*ut.*ux.*(8.*P1.*(ux2+(-1).*uy2)+uy.*((-3).*T3.*(2+3.*ux2).*uy+P2.*ux.*(20+6.*ux2+6.*uy2)));
    Rdtuydt = (1/4).*(8.*P1.*uy.*(utdy.*((-2)+(-4).*uy2)+3.*ut.*uy.*uydy)+ux.*(6.*T2.*Tdx.*ut+4.*P2dy.*ut.*(1+ut2)+16.*P2.*utdy+8.*P1dx.*ut.*(1+ux2)+4.*uy.*((4.*P1+3.*T3).*ut.*uxdy+4.*P2.*(3.*utdy.*uy+ut.*(uxdx+(-1).*uydy)))+ux.*(3.*T3.*ut.*(uxdx+(-1).*uydy)+4.*ut.*((-1).*P2.*(5.*uxdy+uydx)+P1.*(uxdx+uydy))+uy.*(16.*P2dx.*ut+24.*P2.*utdx+(-24).*P1.*utdy+(-18).*T3.*utdy+uy.*((24.*P2.*utdx+(-48).*P1.*utdy+(-18).*T3.*utdy).*uy+27.*T3.*ut.*uydy+(-18).*ut.*(P2.*(uxdy+uydx)+(-2).*P1.*uydy)))+ux.*(12.*P1.*(utdx+utdx.*ux2+(-1).*ut.*ux.*uxdx)+(-9).*T3.*(ut2.*utdx+(-1).*ux.*(ut.*uxdx+(-2).*utdy.*uy)+(-2).*ut.*uy.*(uxdy+uydx))+2.*P2.*(utdy.*(14+3.*ux2+27.*uy2)+3.*(ux.*(6.*utdx.*uy+(-1).*ut.*(uxdy+uydx))+(-4).*ut.*uy.*(uxdx+uydy)))))));
    Rdtuydtdx = (-1/4).*ut.*ux2.*uy.*((-9).*T3.*ux.*uy+P2.*(4+6.*ux2+6.*uy2));
    Rdtuydtdy = (1/4).*ut.*uy.*(8.*P1.*uy2+ux.*((-8).*P2.*uy+ux.*(4.*P1+(-3).*T3+(-12).*P2.*ux.*uy+(12.*P1+9.*T3).*uy2)));
    RdtP1dt = 3.*utdx.*ux3.*(1+ux2).*uy+utdy.*(2.*ux4+ux2.*(2+(-3).*uy2+(-3).*uy4)+(-2).*(uy2+uy4))+ut.*((-2).*ux3.*uxdy+(-3).*ux4.*uxdx.*uy+2.*ux.*uxdy.*uy2+2.*uy3.*uydy+ux2.*uy.*(uxdx+uydy+3.*uy2.*uydy));
    RdtP1dtdx = 2.*ut.*ux.*(1+ux2).*uy;
    RdtP1dtdy = 0;
    RdtP2dt = (1/2).*ux.*uy.*(utdy.*(8+8.*uy2+ux2.*(14+3.*ux2+9.*uy2))+4.*ut.*uy.*(uxdx+(-1).*uydy)+ux.*(3.*utdx.*uy.*(2+3.*ux2+uy2)+ut.*((-10).*uxdy+(-2).*uydx+(-3).*uy2.*(uxdy+uydx)+(-3).*ux.*(ux.*(uxdy+uydx)+2.*uy.*(uxdx+uydy)))));
    RdtP2dtdx = 2.*ut.*ux2.*uy2;
    RdtP2dtdy = ut.*(1+ut2).*ux.*uy;

    Rdt2 = RdtTdt.*Tdt + RdtTdtdx.*Tdtdx + RdtTdtdy.*Tdtdy + Rdtutdt.*utdt + Rdtutdtdx.*utdtdx + Rdtutdtdy.*utdtdy + Rdtuxdt.*uxdt + Rdtuxdtdx.*uxdtdx + Rdtuxdtdy.*uxdtdy + Rdtuydt.*uydt + Rdtuydtdx.*uydtdx + Rdtuydtdy.*uydtdy + RdtP1dt.*P1dt + RdtP1dtdx.*P1dtdx + RdtP1dtdy.*P1dtdy + RdtP2dt.*P2dt + RdtP2dtdx.*P2dtdx + RdtP2dtdy.*P2dtdy;


    RdtTdt = (3/4).*T.*ux.*uy.*(4.*Tdy.*ut+(-3).*T.*(3.*ut2.*utdx.*ux.*uy+3.*ut2.*utdy.*uy2+ut.*(ux.*((-3).*ux.*uxdx.*uy+(-2).*uydx+(-3).*uy2.*(uxdy+uydx))+uy.*(uxdx+((-1)+(-3).*uy2).*uydy))));
    RdtTdtdx = 0;
    RdtTdtdy = (3/2).*T2.*ut.*ux.*uy;
    Rdtutdt = (1/4).*(8.*P1.*uy3.*uydx+ux.*(uy.*((-8).*P1dy+16.*P2dx+6.*T2.*Tdy+uy.*((-4).*P2.*(uxdy+8.*uydx)+2.*uy.*((-4).*P1dy+6.*P2dx+(-3).*P2.*uy.*(uxdy+3.*uydx))+(-4).*P1.*(uxdx+uydy+(-3).*uy2.*uydy)+(-3).*T3.*(uxdx+6.*ut.*utdy.*uy+((-1)+(-3).*uy2).*uydy)))+ux.*(ux.*((-8).*P1.*uxdx+uy.*(12.*P2dx+uy.*((-12).*P1.*uxdx+9.*T3.*uxdx+P2.*((-6).*uxdy+(-18).*uydx))))+uy.*((-8).*P1.*uydx+6.*T3.*uydx+uy.*(8.*P2dy+9.*T3.*((-2).*ut.*utdx+uy.*(uxdy+uydx)))+4.*P2.*((-2).*uxdx+2.*uydy+(-3).*uy2.*(uxdx+uydy))))));
    Rdtutdtdx = (1/4).*(4.*P1.*((-2).*uy2.*(1+uy2)+ux2.*(1+ux2).*(2+3.*uy2))+ux.*uy.*((-9).*T3.*ut2.*ux.*uy+P2.*(16+uy2.*(28+6.*uy2)+ux2.*(16+18.*uy2))));
    Rdtutdtdy = (-3/4).*ux.*uy2.*((-2).*P2.*ux.*(2+ux2+3.*uy2)+uy.*(3.*T3.*ut2+4.*P1.*(1+uy2)));
    Rdtuxdt = (1/4).*(uy.*((-8).*(P1dy+(-1).*P2dx).*ut+6.*T2.*Tdy.*ut+16.*P2.*utdx+uy.*((-4).*P2.*ut.*(uxdy+5.*uydx)+uy.*(((-8).*P1dy+4.*P2dx).*ut+28.*P2.*utdx+(-6).*P2.*uy.*((-1).*utdx.*uy+ut.*(uxdy+uydx)))+(-4).*P1.*(3.*utdy.*uy.*(1+uy2)+ut.*(uxdx+uydy+(-3).*uy2.*uydy))+(-3).*T3.*(3.*utdy.*uy.*(1+uy2)+ut.*(uxdx+((-1)+(-3).*uy2).*uydy))))+ux.*(16.*P1.*utdx+ux.*(4.*P1.*(4.*utdx.*ux+(-3).*ut.*uxdx).*(2+3.*uy2)+uy.*(4.*P2dx.*ut+3.*((-3).*T3.*uy.*(2.*utdx.*ux+(-3).*ut.*uxdx+utdy.*uy)+2.*P2.*(utdx.*(8+9.*uy2)+uy.*(4.*utdy.*ux+(-1).*ut.*(3.*uxdy+uydx))))))+2.*uy.*(6.*T3.*ut.*uydx+(-8).*ut.*(P1.*uydx+P2.*(uxdx+(-1).*uydy))+uy.*(8.*P2dy.*ut+12.*P1.*utdx+(-9).*T3.*utdx+12.*P2.*utdy+3.*uy.*((-3).*(T3.*utdx+(-2).*P2.*utdy).*uy+3.*T3.*ut.*(uxdy+uydx)+(-4).*P2.*ut.*(uxdx+uydy))))));
    Rdtuxdtdx = (1/4).*ut.*ux.*(uy.*(3.*T3.*((-1)+3.*ux2).*uy+P2.*ux.*((-8)+(-12).*uy2))+(-4).*P1.*(uy2+ux2.*(2+3.*uy2)));
    Rdtuxdtdy = (-1/4).*ut.*ux.*uy2.*((-9).*T3.*ux.*uy+P2.*(4+6.*ux2+6.*uy2));
    Rdtuydt = (1/4).*(8.*P1.*uy.*(utdx.*((-2)+(-4).*uy2)+3.*ut.*uy.*uydx)+ux.*(6.*T2.*Tdy.*ut+4.*P2dx.*ut.*(1+ut2)+16.*P2.*utdx+(-8).*P1dy.*ut.*(1+3.*uy2)+uy.*((-6).*T3.*ut.*(uxdx+(-1).*uydy)+(-8).*ut.*(P2.*(uxdy+5.*uydx)+P1.*(uxdx+uydy))+uy.*(84.*P2.*utdx+(-36).*P1.*utdy+(-27).*T3.*utdy+uy.*((30.*P2.*utdx+(-60).*P1.*utdy+(-27).*T3.*utdy).*uy+36.*T3.*ut.*uydy+(-12).*ut.*(P2.*(2.*uxdy+uydx)+(-4).*P1.*uydy))))+ux.*(6.*T3.*ut.*uydx+ux.*(16.*P2.*utdx+(24.*P1.*utdx+(-18).*T3.*utdx+12.*P2.*utdy).*ux.*uy+uy.*(18.*T3.*ut.*uxdx+(54.*P2.*utdx+(-27).*T3.*utdy).*uy+(-12).*ut.*(2.*P1.*uxdx+P2.*(uxdy+uydx))))+(-8).*ut.*(P1.*uydx+P2.*(uxdx+(-1).*uydy))+uy.*(16.*P2dy.*ut+24.*P1.*utdx+(-18).*T3.*utdx+24.*P2.*utdy+9.*uy.*(T3.*((-2).*utdx.*uy+3.*ut.*(uxdy+uydx))+(-4).*P2.*((-2).*utdy.*uy+ut.*(uxdx+uydy)))))));
    Rdtuydtdx = (-1/4).*ut.*uy.*(8.*P1.*(ux2+(-1).*uy2)+ux.*(T3.*ux.*((-6)+(-9).*uy2)+2.*P2.*uy.*(10+3.*ux2+3.*uy2)));
    Rdtuydtdy = (-1/4).*ut.*ux.*uy.*(P2.*ux.*((-8)+12.*uy2)+uy.*(4.*P1+(-3).*T3+((-12).*P1+(-9).*T3).*uy2));
    RdtP1dt = utdx.*((-2).*uy2.*(1+uy2)+ux2.*(1+ux2).*(2+3.*uy2))+2.*ut.*uy3.*uydx+ux.*((-3).*utdy.*uy3.*(1+uy2)+(-1).*ut.*(ux.*(ux.*uxdx.*(2+3.*uy2)+2.*uy.*uydx)+uy2.*(uxdx+uydy+(-3).*uy2.*uydy)));
    RdtP1dtdx = 0;
    RdtP1dtdy = (-2).*ut.*ux.*uy.*(1+uy2);
    RdtP2dt = (1/2).*ux.*uy.*(utdx.*(8+14.*uy2+3.*uy4+ux2.*(8+9.*uy2))+ut.*uy.*((-2).*uxdy+(-10).*uydx+(-3).*uy2.*(uxdy+uydx))+ux.*(3.*utdy.*uy.*(2+ux2+3.*uy2)+ut.*((-4).*uxdx+(-3).*ux.*uy.*(uxdy+uydx)+4.*uydy+(-6).*uy2.*(uxdx+uydy))));
    RdtP2dtdx = ut.*(1+ut2).*ux.*uy;
    RdtP2dtdy = 2.*ut.*ux2.*uy2;

    Rdt3 = RdtTdt.*Tdt + RdtTdtdx.*Tdtdx + RdtTdtdy.*Tdtdy + Rdtutdt.*utdt + Rdtutdtdx.*utdtdx + Rdtutdtdy.*utdtdy + Rdtuxdt.*uxdt + Rdtuxdtdx.*uxdtdx + Rdtuxdtdy.*uxdtdy + Rdtuydt.*uydt + Rdtuydtdx.*uydtdx + Rdtuydtdy.*uydtdy + RdtP1dt.*P1dt + RdtP1dtdx.*P1dtdx + RdtP1dtdy.*P1dtdy + RdtP2dt.*P2dt + RdtP2dtdx.*P2dtdx + RdtP2dtdy.*P2dtdy;


    RdtTdt = (1/4).*T2.*(2.*P1.*(16.*T.*ut+8.*T.*ut.*ux2+8.*T.*ut.*uy2+(-3).*A2.*ux.*uy.*((-1).*ux2+uy2).*((-1).*utdy.*ux+utdx.*uy+ut.*(uxdy+(-1).*uydx))+3.*A1.*(2.*utdy.*uy.*(2+3.*ux4+5.*uy2+3.*uy4+ux2.*(5+2.*uy2))+utdx.*ux.*(4+ux2.*(10+6.*ux2+4.*uy2)+uy2.*(10+6.*uy2))+(-1).*ut.*((-6).*uxdx+(-6).*uydy+uy2.*((-3).*uxdx+uydy+6.*uy2.*uydy)+ux.*(2.*uy.*(uxdy.*(2+uy2)+(2+3.*uy2).*uydx)+ux.*(uxdx+2.*ux.*(3.*ux.*uxdx+uy.*(3.*uxdy+uydx))+(-3).*uydy+2.*uy2.*(uxdx+uydy))))))+3.*(A2.*P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+4.*T3.*ut.*(uxdx+(-1).*uydy)+2.*T3.*uy.*(utdy.*(2+3.*uy2+uy4)+ut.*uy.*(uxdx+(-1).*(3+uy2).*uydy))+(-2).*T3.*ux.*(utdx.*(2+3.*ux2+ux4+uy2+(-1).*uy4)+ut.*uy3.*(uxdy+uydx)+ux.*(utdy.*((-1)+ux2).*uy+(-1).*ut.*(uxdx.*(3+ux2+(-1).*uy2)+ux.*uy.*(uxdy+uydx)+((-1)+uy2).*uydy)))+4.*A1.*(ut.*uy.*(P1dy+P2.*uxdy.*uy).*(2+uy2)+ux.*(P1dx.*ut.*(1+ut2)+P2.*uy.*((-2).*utdy.*uy3+ut.*(uxdx.*(2+uy2)+((-2)+uy2).*uydy))+ux.*(P1dy.*ut.*uy+P2.*(2.*utdy.*ux.*uy2+(-2).*utdx.*uy.*((-1).*ux2+uy2)+(-1).*ut.*(uy2.*(uxdy+(-1).*uydx)+2.*uydx+ux.*(ux.*uydx+uy.*(uxdx+uydy)))))))));
    RdtTdtdx = 0;
    RdtTdtdy = 0;
    Rdtutdt = (1/12).*((-16).*A1.*P1.^2.*(2.*ut.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2))+(-1).*(4+3.*ux2+3.*uy2).*((1+ux2).*uxdx+(-1).*(1+uy2).*uydy))+2.*P1.*(T4.*(24+18.*ux2+18.*uy2)+4.*A1.*P2.*(12.*uxdy+(-2).*ut.*(utdx.*uy.*(2+3.*ux2+uy2)+utdy.*ux.*(2+ux2+3.*uy2))+12.*uydx+uy2.*(16+5.*uy2).*(uxdy+uydx)+ux.*(5.*ux3.*(uxdy+uydx)+2.*ux.*(8+5.*uy2).*(uxdy+uydx)+6.*ux2.*uy.*(uxdx+uydy)+2.*uy.*(4+3.*uy2).*(uxdx+uydy)))+3.*T3.*(A2.*ux.*uy.*(ux2+(-1).*uy2).*(uxdy+(-1).*uydx)+(-1).*A1.*((-6).*uxdx+(-6).*uydy+uy2.*((-3).*uxdx+uydy+6.*uy2.*uydy)+ux.*(2.*uy.*(uxdy.*(2+uy2)+(2+3.*uy2).*uydx)+ux.*(uxdx+2.*ux.*(3.*ux.*uxdx+uy.*(3.*uxdy+uydx))+(-3).*uydy+2.*uy2.*(uxdx+uydy))))))+3.*T3.*((-1).*A2.*P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2).*(uxdy+(-1).*uydx)+2.*T3.*(uxdx+(-1).*uydy)+T3.*uy2.*(uxdx+(-1).*(3+uy2).*uydy)+T3.*ux.*((-1).*uy3.*(uxdy+uydx)+ux.*(uxdx.*(3+ux2+(-1).*uy2)+ux.*uy.*(uxdy+uydx)+((-1)+uy2).*uydy))+4.*A1.*(P1dx.*ux.*(4+3.*ux2+3.*uy2)+P1dy.*uy.*(4+3.*ux2+3.*uy2)+P2.*(uxdy.*uy2.*(2+uy2)+(-1).*ux.*((-1).*uy.*(uxdx.*(2+uy2)+((-2)+uy2).*uydy)+ux.*(uy2.*(uxdy+(-1).*uydx)+2.*uydx+ux.*(ux.*uydx+uy.*(uxdx+uydy))))))));
    Rdtutdtdx = (1/12).*((-3).*T.^6.*ux.*(2+3.*ux2+ux4+uy2+(-1).*uy4)+(-3).*A2.*T3.*uy.*((-2).*P1.*ux.*uy.*(ux2+(-1).*uy2)+P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2))+4.*A1.*((-4).*P1.^2.*(1+ut2).*ux.*(1+ux2)+6.*P2.*T3.*ux2.*uy.*(ux2+(-1).*uy2)+P1.*((-2).*P2.*uy.*(3.*ux4+4.*ux2.*(2+uy2)+(2+uy2).^2)+T3.*ux.*(6+ux2.*(15+9.*ux2+6.*uy2)+uy2.*(15+9.*uy2)))));
    Rdtutdtdy = (1/12).*(3.*T3.*(T3.*uy.*(2+ux2+(-1).*ux4+3.*uy2+uy4)+A2.*ux.*((-2).*P1.*ux.*uy.*(ux2+(-1).*uy2)+P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2)))+4.*A1.*(6.*P2.*T3.*ux.*uy2.*(ux2+(-1).*uy2)+4.*P1.^2.*(1+ut2).*uy.*(1+uy2)+P1.*((-2).*P2.*ux.*(4+ux4+8.*uy2+3.*uy4+4.*ux2.*(1+uy2))+3.*T3.*uy.*(2+3.*ux4+5.*uy2+3.*uy4+ux2.*(5+2.*uy2)))));
    Rdtuxdt = (1/12).*((-6).*T.^6.*utdx+3.*A2.*T3.*(P2.*(utdy.*(5.*ux4+(-6).*ux2.*((-2)+uy2)+(2+uy2).^2)+(-4).*ux.*(2+ux2+(-1).*uy2).*(utdx.*uy+ut.*(uxdy+(-1).*uydx)))+(-2).*P1.*uy.*(utdy.*ux.*(4.*ux2+(-2).*uy2)+((-3).*ux2+uy2).*(utdx.*uy+ut.*(uxdy+(-1).*uydx))))+3.*T.^6.*uy2.*(utdx.*((-1)+uy2)+(-1).*ut.*uy.*(uxdy+uydx))+(-3).*T.^6.*ux.*(utdx.*ux.*(9+5.*ux2)+2.*utdy.*((-1)+2.*ux2).*uy+ut.*(uxdx.*((-6)+(-4).*ux2+2.*uy2)+(-3).*ux.*uy.*(uxdy+uydx)+(-2).*((-1)+uy2).*uydy))+4.*A1.*(3.*P1dx.*T3.*ut.*(1+ut2)+(-4).*P1.^2.*(1+ut2).*(utdx+3.*utdx.*ux2+(-2).*ut.*ux.*uxdx)+3.*P2.*T3.*((-2).*utdy.*uy2.*((-3).*ux2+uy2)+(-4).*utdx.*ux.*uy.*((-2).*ux2+uy2)+ut.*(uy.*(uxdx.*(2+uy2)+((-2)+uy2).*uydy)+ux.*((-2).*uy2.*(uxdy+(-1).*uydx)+(-4).*uydx+ux.*((-4).*ux.*uydx+(-3).*uy.*(uxdx+uydy)))))+P1.*((-2).*P2.*(1+ut2).*(utdy.*(2+3.*ux2+3.*uy2)+2.*uy.*(3.*utdx.*ux+(-1).*ut.*(uxdx+uydy)))+T3.*(utdx.*(6+uy2.*(15+9.*uy2)+ux2.*(45+45.*ux2+18.*uy2))+3.*(2.*utdy.*ux.*uy.*(5+6.*ux2+2.*uy2)+(-1).*ut.*(uy.*(uxdy.*(2+uy2)+(2+3.*uy2).*uydx)+ux.*(uxdx+3.*ux.*(4.*ux.*uxdx+uy.*(3.*uxdy+uydx))+(-3).*uydy+2.*uy2.*(uxdx+uydy))))))));
    Rdtuxdtdx = (1/12).*ut.*(3.*T.^6.*(2+ux4+uy2+(-1).*ux2.*((-3)+uy2))+2.*A1.*(8.*P1.^2.*(1+ut2).*(1+ux2)+(-6).*P2.*T3.*ux.*uy.*((-2)+ux2+(-1).*uy2)+P1.*(8.*P2.*(1+ut2).*ux.*uy+(-3).*T3.*((-6)+ux2+6.*ux4+(-3).*uy2+2.*ux2.*uy2))));
    Rdtuxdtdy = (1/12).*ut.*(4.*A1.*(3.*P2.*T3.*uy2.*(2+(-1).*ux2+uy2)+P1.*(2.*P2.*(1+ut2).^2+(-3).*T3.*ux.*uy.*(2+3.*ux2+uy2)))+(-3).*T3.*(T3.*ux.*uy.*((-1).*ux2+uy2)+A2.*((-2).*P1.*ux.*uy.*(ux2+(-1).*uy2)+P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2))));
    Rdtuydt = (1/12).*(6.*T.^6.*utdy+3.*A2.*T3.*((-2).*P1.*ux.*(utdy.*ux.*(ux2+(-3).*uy2)+(-2).*utdx.*uy.*(ux2+(-2).*uy2)+(-1).*ut.*(ux2+(-3).*uy2).*(uxdy+(-1).*uydx))+(-1).*P2.*(utdx.*(4+ux4+12.*uy2+5.*uy4+ux2.*(4+(-6).*uy2))+(-4).*uy.*(2+(-1).*ux2+uy2).*(utdy.*ux+ut.*((-1).*uxdy+uydx))))+3.*T.^6.*ux.*(uy.*(utdx.*((-2)+4.*uy2)+(-3).*ut.*uy.*(uxdy+uydx))+ux.*((-1).*utdy.*((-1)+ux2)+ut.*(ux.*(uxdy+uydx)+(-2).*uy.*(uxdx+(-1).*uydy))))+3.*T.^6.*uy.*(utdy.*uy.*(9+5.*uy2)+2.*ut.*(uxdx+((-3)+(-2).*uy2).*uydy))+A1.*(12.*P1dy.*T3.*ut.*(1+ut2)+24.*P1.*T3.*utdy+32.*P1.*((-1).*P2.*utdx+P1.*utdy)+uy.*(48.*P2.*T3.*ut.*uxdy.*(1+uy2)+(-16).*P1.^2.*(utdy.*uy.*((-7)+(-3).*uy2)+2.*ut.*(2+uy2).*uydy)+P1.*((-8).*P2.*utdx.*uy.*(8+3.*uy2)+T3.*(180.*utdy.*uy.*(1+uy2)+ut.*(36.*uxdx+(-12).*uydy+(-144).*uy2.*uydy))))+4.*ux.*(4.*P1.^2.*ux.*(utdy+3.*utdy.*uy2+(-2).*ut.*uy.*uydy)+3.*P2.*T3.*(2.*utdx.*ux.*(ux2+(-3).*uy2)+4.*utdy.*uy.*(ux2+(-2).*uy2)+(-1).*ut.*((-2).*uxdx+2.*uydy+(-3).*uy2.*(uxdx+uydy)+ux.*(2.*uy.*(uxdy+(-1).*uydx)+ux.*(uxdx+uydy))))+P1.*(2.*P2.*(utdx.*ux.*((-8)+(-3).*ux2+(-6).*uy2)+2.*(1+ut2).*((-3).*utdy.*uy+ut.*(uxdx+uydy)))+3.*T3.*(2.*utdx.*uy.*(5+2.*ux2+6.*uy2)+utdy.*ux.*(5+3.*ux2+6.*uy2)+(-1).*ut.*(2.*uxdy+2.*uydx+3.*uy2.*(uxdy+3.*uydx)+ux.*(ux.*(3.*uxdy+uydx)+2.*uy.*(uxdx+uydy))))))));
    Rdtuydtdx = (1/12).*ut.*(4.*A1.*((-3).*P2.*T3.*ux2.*(2+ux2+(-1).*uy2)+P1.*(2.*P2.*(1+ut2).^2+(-3).*T3.*ux.*uy.*(2+ux2+3.*uy2)))+3.*T3.*(T3.*ux.*uy.*(ux2+(-1).*uy2)+A2.*((-2).*P1.*ux.*uy.*(ux2+(-1).*uy2)+P2.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2))));
    Rdtuydtdy = (1/12).*ut.*((-3).*T.^6.*(2+3.*uy2+uy4+(-1).*ux2.*((-1)+uy2))+(-2).*A1.*(6.*P2.*T3.*ux.*uy.*(2+ux2+(-1).*uy2)+8.*P1.^2.*(1+ut2).*(1+uy2)+P1.*((-8).*P2.*(1+ut2).*ux.*uy+3.*T3.*((-6)+uy2+6.*uy4+ux2.*((-3)+2.*uy2)))));
    RdtP1dt = (1/6).*(6.*T4.*ut.*(1+ut2)+(-4).*A1.*(1+ut2).*(4.*P1.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy))+P2.*(utdx.*uy.*(2+3.*ux2+uy2)+utdy.*ux.*(2+ux2+3.*uy2)+(-1).*ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+2.*ux.*uy.*(uxdx+uydy))))+3.*T3.*((-1).*A2.*ux.*uy.*(ux2+(-1).*uy2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+A1.*(2.*utdy.*uy.*(2+3.*ux4+5.*uy2+3.*uy4+ux2.*(5+2.*uy2))+utdx.*ux.*(4+ux2.*(10+6.*ux2+4.*uy2)+uy2.*(10+6.*uy2))+(-1).*ut.*((-6).*uxdx+(-6).*uydy+uy2.*((-3).*uxdx+uydy+6.*uy2.*uydy)+ux.*(2.*uy.*(uxdy.*(2+uy2)+(2+3.*uy2).*uydx)+ux.*(uxdx+2.*ux.*(3.*ux.*uxdx+uy.*(3.*uxdy+uydx))+(-3).*uydy+2.*uy2.*(uxdx+uydy)))))));
    RdtP1dtdx = A1.*T3.*ut.*(1+ut2).*ux;
    RdtP1dtdy = A1.*T3.*ut.*(1+ut2).*uy;
    RdtP2dt = (1/12).*(3.*A2.*T3.*(ux4+(-2).*ux2.*((-2)+uy2)+(2+uy2).^2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+4.*A1.*(2.*P1.*(1+ut2).*((-1).*utdx.*uy.*(2+3.*ux2+uy2)+(-1).*utdy.*ux.*(2+ux2+3.*uy2)+ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+2.*ux.*uy.*(uxdx+uydy)))+3.*T3.*(ut.*uxdy.*uy2.*(2+uy2)+ux.*(uy.*((-2).*utdy.*uy3+ut.*(uxdx.*(2+uy2)+((-2)+uy2).*uydy))+ux.*(2.*utdy.*ux.*uy2+(-2).*utdx.*uy.*((-1).*ux2+uy2)+(-1).*ut.*(uy2.*(uxdy+(-1).*uydx)+2.*uydx+ux.*(ux.*uydx+uy.*(uxdx+uydy))))))));
    RdtP2dtdx = 0;
    RdtP2dtdy = 0;

    Rdt4 = RdtTdt.*Tdt + RdtTdtdx.*Tdtdx + RdtTdtdy.*Tdtdy + Rdtutdt.*utdt + Rdtutdtdx.*utdtdx + Rdtutdtdy.*utdtdy + Rdtuxdt.*uxdt + Rdtuxdtdx.*uxdtdx + Rdtuxdtdy.*uxdtdy + Rdtuydt.*uydt + Rdtuydtdx.*uydtdx + Rdtuydtdy.*uydtdy + RdtP1dt.*P1dt + RdtP1dtdx.*P1dtdx + RdtP1dtdy.*P1dtdy + RdtP2dt.*P2dt + RdtP2dtdx.*P2dtdx + RdtP2dtdy.*P2dtdy;


    RdtTdt = (1/2).*T2.*(ux.*uy.*(P2.*(8.*T.*ut.*(1+ut2)+3.*A2.*ux.*uy.*(ux2+(-1).*uy2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx)))+6.*((-1).*A2.*P1.*(1+ux2).*(1+uy2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+(-1).*T3.*(ut2.*utdx.*(1+ux2).*uy+ut2.*utdy.*ux.*(1+uy2)+(-1).*ut.*(uxdy+uydx+uy2.*(uxdy+uydx)+ux.*(ux.*(uxdy+ux.*uxdx.*uy+uydx+uy2.*(uxdy+uydx))+uy.*(uxdx+uydy+uy2.*uydy))))))+3.*A1.*(4.*P1.*(utdx.*ux.*(1+ux2).*(ux2+(-1).*uy2).*(1+uy2)+utdy.*(1+ux2).*uy.*(ux2+(-1).*uy2).*(1+uy2)+ut.*(uy4.*uydy+ux.*(uy3.*uydx+ux.*((-1).*ux.*(ux.*uxdx.*(1+uy2)+uy.*(uxdy+uy2.*(uxdy+(-1).*uydx)))+uy4.*uydy))))+ux.*uy.*(2.*P2dx.*ut.*(1+ut2).*ux+2.*P2dy.*ut.*(1+ut2).*uy+P2.*(2.*utdy.*uy.*(6+2.*ux4+9.*uy2+2.*uy4+ux2.*(9+8.*uy2))+utdx.*ux.*(12+uy2.*(18+4.*uy2)+ux2.*(18+4.*ux2+16.*uy2))+ut.*(6.*(uxdx+uydy)+uy2.*(3.*uxdx+((-9)+(-4).*uy2).*uydy)+ux.*((-4).*uy.*(uxdy.*(3+2.*uy2)+(3+uy2).*uydx)+(-1).*ux.*(9.*uxdx+4.*ux.*(ux.*uxdx+uy.*(uxdy+2.*uydx))+(-3).*uydy+8.*uy2.*(uxdx+uydy))))))));
    RdtTdtdx = 0;
    RdtTdtdy = 0;
    Rdtutdt = (1/6).*(3.*A2.*T3.*ux.*uy.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)).*(uxdy+(-1).*uydx)+3.*T4.*ux.*uy.*(P2.*(8+6.*ux2+6.*uy2)+T2.*(uxdy+(-2).*ut.*(utdx.*(1+ux2).*uy+utdy.*ux.*(1+uy2))+uydx+uy2.*(uxdy+uydx)+ux.*(ux.*(uxdy+ux.*uxdx.*uy+uydx+uy2.*(uxdy+uydx))+uy.*(uxdx+uydy+uy2.*uydy))))+A1.*(12.*P1.*T3.*uy4.*uydy+ux.*(uy.*(6.*P2dy.*T3.*uy.*(4+3.*uy2)+4.*P2.^2.*(2+uy2).*(6.*uxdy+6.*uydx+uy.*((-2).*ut.*utdx+5.*uy.*(uxdy+uydx)))+3.*P2.*T3.*(3.*uxdx.*(2+uy2)+(6+(-9).*uy2+(-4).*uy4).*uydy)+P1.*(12.*T3.*uy2.*uydx+P2.*(uxdx.*(32+24.*uy2)+(-8).*(1+uy2).*((-2).*ut.*utdy.*uy+(4+3.*uy2).*uydy))))+ux.*(ux.*(4.*P1.*((-3).*T3.*(ux.*uxdx.*(1+uy2)+uy.*(uxdy+uy2.*(uxdy+(-1).*uydx)))+2.*P2.*uy.*((-2).*ut.*utdx.*ux+uxdx.*(7+3.*ux2+3.*uy2)+(-3).*(1+uy2).*uydy))+uy.*(18.*P2dx.*T3.*ux+18.*P2dy.*T3.*uy+(-3).*P2.*T3.*(9.*uxdx+4.*ux.*(ux.*uxdx+uy.*(uxdy+2.*uydx))+(-3).*uydy+8.*uy2.*(uxdx+uydy))+4.*P2.^2.*(16.*uxdy+(-2).*ut.*(utdy.*ux+3.*utdx.*uy)+16.*uydx+10.*uy2.*(uxdy+uydx)+ux.*(5.*ux.*(uxdy+uydx)+6.*uy.*(uxdx+uydy)))))+2.*uy.*(12.*P2dx.*T3+(-8).*P2.*ut.*(P1.*utdx+P2.*utdy)+uy.*((-18).*P2.*T3.*(uxdy+uydx)+16.*P2.^2.*(uxdx+uydy)+3.*uy.*(3.*P2dx.*T3+(-4).*P2.^2.*ut.*utdy+2.*uy.*(2.*P2.^2.*(uxdx+uydy)+(-1).*T3.*(P2.*(2.*uxdy+uydx)+(-1).*P1.*uydy)))))))));
    Rdtutdtdx = (1/6).*ux.*((-3).*T3.*uy2.*(T3.*ut2.*(1+ux2)+A2.*(P2.*ux.*uy.*(ux2+(-1).*uy2)+(-2).*P1.*(1+ux2).*(1+uy2)))+2.*A1.*((-2).*P1.*(1+ux2).*(2.*P2.*(1+ut2).*ux.*uy+(-3).*T3.*(ux2+(-1).*uy2).*(1+uy2))+P2.*uy.*((-2).*P2.*uy.*(3.*ux4+4.*ux2.*(2+uy2)+(2+uy2).^2)+T3.*ux.*(18+uy2.*(27+6.*uy2)+ux2.*(27+6.*ux2+24.*uy2)))));
    Rdtutdtdy = (1/6).*uy.*((-3).*T3.*ux2.*(T3.*ut2.*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2)))+(-2).*A1.*((-6).*P1.*T3.*(1+ux2).*(ux2+(-1).*uy2).*(1+uy2)+2.*P2.^2.*ux2.*(4+ux4+8.*uy2+3.*uy4+4.*ux2.*(1+uy2))+(-1).*P2.*ux.*uy.*(4.*P1.*(1+ut2).*(1+uy2)+3.*T3.*(6+2.*ux4+9.*uy2+2.*uy4+ux2.*(9+8.*uy2)))));
    Rdtuxdt = (1/6).*(3.*T3.*uy.*(P2.*(2.*T.*ut.*(1+ut2)+A2.*ux.*uy.*(utdy.*ux.*(5.*ux2+(-3).*uy2)+2.*((-2).*ux2+uy2).*(utdx.*uy+ut.*(uxdy+(-1).*uydx))))+2.*A2.*P1.*(1+uy2).*((-2).*utdy.*(ux+2.*ux3)+(1+3.*ux2).*(utdx.*uy+ut.*(uxdy+(-1).*uydx)))+(-1).*T3.*(ut2.*utdx.*(1+3.*ux2).*uy+2.*ut2.*utdy.*ux.*(1+uy2)+(-1).*ut.*(uxdy+uydx+uy2.*(uxdy+uydx)+ux.*(4.*ux2.*uxdx.*uy+3.*ux.*(1+uy2).*(uxdy+uydx)+2.*uy.*(uxdx+uydy+uy2.*uydy)))))+A1.*(4.*P1.*(2.*P2.*(1+ut2).*uy.*((-2).*utdx.*(ux+2.*ux3)+utdy.*uy.*(1+uy2)+ut.*(uxdx+3.*ux2.*uxdx+(-1).*(1+uy2).*uydy))+3.*T3.*(utdx.*(1+uy2).*(5.*ux4+(-1).*uy2+(-3).*ux2.*((-1)+uy2))+ut.*uy3.*uydx+ux.*((-2).*utdy.*uy.*(1+uy2).*((-1)+(-2).*ux2+uy2)+ut.*(ux.*((-4).*ux.*uxdx.*(1+uy2)+(-3).*uy.*(uxdy+uy2.*(uxdy+(-1).*uydx)))+2.*uy4.*uydy))))+uy.*(12.*P2dx.*T3.*ut.*(1+ut2).*ux+6.*P2dy.*T3.*ut.*(1+ut2).*uy+(-4).*P2.^2.*(1+ut2).*(utdx.*uy.*(2+9.*ux2+uy2)+utdy.*ux.*(4+4.*ux2+6.*uy2)+(-1).*ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+4.*ux.*uy.*(uxdx+uydy)))+3.*P2.*T3.*(2.*utdy.*uy.*(6+10.*ux4+9.*uy2+2.*uy4+ux2.*(27+24.*uy2))+utdx.*ux.*(24+uy2.*(36+8.*uy2)+ux2.*(72+24.*ux2+64.*uy2))+ut.*(6.*(uxdx+uydy)+uy2.*(3.*uxdx+((-9)+(-4).*uy2).*uydy)+ux.*((-20).*ux3.*uxdx+(-16).*ux2.*uy.*(uxdy+2.*uydx)+(-8).*uy.*(uxdy.*(3+2.*uy2)+(3+uy2).*uydx)+3.*ux.*((-9).*uxdx+3.*uydy+(-8).*uy2.*(uxdx+uydy))))))));
    Rdtuxdtdx = (1/6).*ut.*ux.*(3.*T.^6.*ux.*(1+ux2).*uy2+A1.*(4.*P1.*(2.*P2.*(1+ut2).*(1+ux2).*uy+(-3).*T3.*ux3.*(1+uy2))+P2.*uy.*(8.*P2.*(1+ut2).*ux.*uy+T3.*(18+9.*uy2+ux2.*((-27)+(-12).*ux2+(-24).*uy2)))));
    Rdtuxdtdy = (1/6).*ut.*ux.*uy.*(4.*A1.*(P2.^2.*(1+ut2).^2+(-3).*P1.*T3.*ux2.*(1+uy2)+(-3).*P2.*T3.*ux.*uy.*(3+ux2+2.*uy2))+3.*T3.*(T3.*(1+ux2).*(1+uy2)+A2.*(2.*P1.*(1+ux2).*(1+uy2)+P2.*ux.*uy.*((-1).*ux2+uy2))));
    Rdtuydt = (1/6).*(4.*A1.*P2.^2.*(1+ut2).*ux.*((-2).*utdx.*uy.*(2+3.*ux2+2.*uy2)+(-1).*utdy.*ux.*(2+ux2+9.*uy2)+ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+4.*ux.*uy.*(uxdx+uydy)))+P2.*ux.*(6.*T4.*ut.*(1+ut2)+(-8).*A1.*P1.*(1+ut2).*(utdx.*(ux+ux3)+(-2).*utdy.*uy.*(1+2.*uy2)+(-1).*ut.*((1+ux2).*uxdx+((-1)+(-3).*uy2).*uydy))+3.*T3.*(A2.*ux.*uy.*(2.*utdy.*ux.*(ux2+(-2).*uy2)+utdx.*uy.*((-3).*ux2+5.*uy2)+(-2).*ut.*(ux2+(-2).*uy2).*(uxdy+(-1).*uydx))+A1.*(4.*utdy.*uy.*(6+6.*uy2.*(3+uy2)+ux2.*(9+2.*ux2+16.*uy2))+utdx.*ux.*(12+uy2.*(54+20.*uy2)+ux2.*(18+4.*ux2+48.*uy2))+ut.*(6.*(uxdx+uydy)+uy2.*(9.*uxdx+((-27)+(-20).*uy2).*uydy)+ux.*((-4).*ux3.*uxdx+(-8).*ux2.*uy.*(uxdy+2.*uydx)+8.*uy.*(uxdy.*((-3)+(-4).*uy2)+((-3)+(-2).*uy2).*uydx)+3.*ux.*((-3).*uxdx+uydy+(-8).*uy2.*(uxdx+uydy)))))))+3.*T3.*((-2).*A2.*P1.*ux.*(1+ux2).*((-2).*utdx.*uy.*(1+2.*uy2)+utdy.*ux.*(1+3.*uy2)+(-1).*ut.*(1+3.*uy2).*(uxdy+(-1).*uydx))+T3.*ux.*((-2).*ut2.*utdx.*(1+ux2).*uy+(-1).*ut2.*utdy.*ux.*(1+3.*uy2)+ut.*(uxdy+uydx+3.*uy2.*(uxdy+uydx)+ux.*(ux.*(uxdy+2.*ux.*uxdx.*uy+uydx+3.*uy2.*(uxdy+uydx))+2.*uy.*(uxdx+uydy+2.*uy2.*uydy))))+2.*A1.*(2.*P1.*uy2.*(utdy.*((-3)+(-5).*uy2)+4.*ut.*uy.*uydy)+ux.*(P2dx.*ut.*(1+ut2).*ux+2.*(P2dy.*ut.*(1+ut2).*uy+P1.*(2.*utdx.*(1+ux2).*uy.*((-1)+ux2+(-2).*uy2)+utdy.*ux.*(1+ux2+3.*ux2.*uy2+(-5).*uy4)+(-1).*ut.*((-3).*uy2.*uydx+ux.*(ux.*(uxdy+2.*ux.*uxdx.*uy+3.*uy2.*(uxdy+(-1).*uydx))+(-4).*uy3.*uydy))))))));
    Rdtuydtdx = (1/6).*ut.*ux.*uy.*(4.*A1.*(P2.^2.*(1+ut2).^2+3.*P1.*T3.*(1+ux2).*uy2+(-3).*P2.*T3.*ux.*uy.*(3+2.*ux2+uy2))+3.*T3.*(T3.*(1+ux2).*(1+uy2)+A2.*(P2.*ux.*uy.*(ux2+(-1).*uy2)+(-2).*P1.*(1+ux2).*(1+uy2))));
    Rdtuydtdy = (1/6).*ut.*uy.*(3.*T.^6.*ux2.*uy.*(1+uy2)+A1.*(4.*P1.*(3.*T3.*(1+ux2).*uy3+(-2).*P2.*(1+ut2).*ux.*(1+uy2))+P2.*ux.*(8.*P2.*(1+ut2).*ux.*uy+T3.*(18+ux2.*(9+(-24).*uy2)+uy2.*((-27)+(-12).*uy2)))));
    RdtP1dt = (1/3).*((-3).*A2.*T3.*ux.*(1+ux2).*uy.*(1+uy2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+A1.*(6.*T3.*(utdx.*(ux+ux3).*(ux2+(-1).*uy2).*(1+uy2)+utdy.*(1+ux2).*uy.*(ux2+(-1).*uy2).*(1+uy2)+ut.*((-1).*ux4.*uxdx.*(1+uy2)+ux.*uy3.*uydx+ux3.*uy.*((-1).*uxdy+uy2.*((-1).*uxdy+uydx))+uy4.*uydy+ux2.*uy4.*uydy))+(-4).*P2.*(1+ut2).*ux.*uy.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy))));
    RdtP1dtdx = 0;
    RdtP1dtdy = 0;
    RdtP2dt = (1/6).*ux.*uy.*(6.*T4.*ut.*(1+ut2)+(-8).*A1.*(1+ut2).*(P1.*(utdx.*(ux+ux3)+(-1).*utdy.*uy.*(1+uy2)+ut.*((-1).*(1+ux2).*uxdx+(1+uy2).*uydy))+P2.*(utdx.*uy.*(2+3.*ux2+uy2)+utdy.*ux.*(2+ux2+3.*uy2)+(-1).*ut.*((1+ut2).*uxdy+ux2.*uydx+(2+uy2).*uydx+2.*ux.*uy.*(uxdx+uydy))))+3.*T3.*(A2.*ux.*uy.*(ux2+(-1).*uy2).*(utdy.*ux+(-1).*utdx.*uy+ut.*((-1).*uxdy+uydx))+A1.*(2.*utdy.*uy.*(6+2.*ux4+9.*uy2+2.*uy4+ux2.*(9+8.*uy2))+utdx.*ux.*(12+uy2.*(18+4.*uy2)+ux2.*(18+4.*ux2+16.*uy2))+ut.*(6.*(uxdx+uydy)+uy2.*(3.*uxdx+((-9)+(-4).*uy2).*uydy)+ux.*((-4).*uy.*(uxdy.*(3+2.*uy2)+(3+uy2).*uydx)+(-1).*ux.*(9.*uxdx+4.*ux.*(ux.*uxdx+uy.*(uxdy+2.*uydx))+(-3).*uydy+8.*uy2.*(uxdx+uydy)))))));
    RdtP2dtdx = A1.*T3.*ut.*(1+ut2).*ux2.*uy;
    RdtP2dtdy = A1.*T3.*ut.*(1+ut2).*ux.*uy2;

    Rdt5 = RdtTdt.*Tdt + RdtTdtdx.*Tdtdx + RdtTdtdy.*Tdtdy + Rdtutdt.*utdt + Rdtutdtdx.*utdtdx + Rdtutdtdy.*utdtdy + Rdtuxdt.*uxdt + Rdtuxdtdx.*uxdtdx + Rdtuxdtdy.*uxdtdy + Rdtuydt.*uydt + Rdtuydtdx.*uydtdx + Rdtuydtdy.*uydtdy + RdtP1dt.*P1dt + RdtP1dtdx.*P1dtdx + RdtP1dtdy.*P1dtdy + RdtP2dt.*P2dt + RdtP2dtdx.*P2dtdx + RdtP2dtdy.*P2dtdy;


    xt.toc('rp');

    RdtTdt    = []; 
    RdtTdtdx  = []; 
    RdtTdtdy  = []; 
    Rdtutdt   = []; 
    Rdtutdtdx = []; 
    Rdtutdtdy = []; 
    Rdtuxdt   = []; 
    Rdtuxdtdx = []; 
    Rdtuxdtdy = []; 
    Rdtuydt   = []; 
    Rdtuydtdx = []; 
    Rdtuydtdy = []; 
    RdtP1dt   = []; 
    RdtP1dtdx = []; 
    RdtP1dtdy = []; 
    RdtP2dt   = []; 
    RdtP2dtdx = []; 
    RdtP2dtdy = []; 



    R1 = Rdt1 + (mdt11 .* Tdt + mdt12 .* uxdt + mdt13 .* uydt);
    R2 = Rdt2 + (               mdt22 .* uxdt + mdt23 .* uydt + mdt24 .* P1dt + mdt25 .* P2dt);
    R3 = Rdt3 + (               mdt32 .* uxdt + mdt33 .* uydt + mdt34 .* P1dt + mdt35 .* P2dt);
    R4 = Rdt4 + (               mdt42 .* uxdt + mdt43 .* uydt + mdt44 .* P1dt                );
    R5 = Rdt5 + (               mdt52 .* uxdt + mdt53 .* uydt +                 mdt55 .* P2dt);





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

    
    Tdtdt  = -(i11.*R1 + i12.*R2 + i13.*R3 + i14.*R4 + i15.*R5) ./ mdet;
    uxdtdt = -(          i22.*R2 + i23.*R3 + i24.*R4 + i25.*R5) ./ mdet;
    uydtdt = -(          i32.*R2 + i33.*R3 + i34.*R4 + i35.*R5) ./ mdet;
    P1dtdt = -(          i42.*R2 + i43.*R3 + i44.*R4 + i45.*R5) ./ mdet;
    P2dtdt = -(          i52.*R2 + i53.*R3 + i54.*R4 + i55.*R5) ./ mdet;

    if ...
      any(~isfinite( Tdtdt)) | ...
      any(~isfinite(uxdtdt)) | ...
      any(~isfinite(uydtdt)) | ...
      any(~isfinite(P1dtdt)) | ...
      any(~isfinite(P2dtdt)) 
      error('ecountered infinity, got scared, and run away screaming.');
    end





    udtdt = zeros(Nx.*Ny.*5, 1);

    udtdt(0*Nt+1 : 1*Nt) =  Tdtdt(:);
    udtdt(1*Nt+1 : 2*Nt) = uxdtdt(:);
    udtdt(2*Nt+1 : 3*Nt) = uydtdt(:);
    udtdt(3*Nt+1 : 4*Nt) = P1dtdt(:);
    udtdt(4*Nt+1 : 5*Nt) = P2dtdt(:);


    

end
