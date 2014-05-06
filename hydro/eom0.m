function udt = eom0(t, u, params)
%   choice of Pi: 
%   Pi_xx =  2 P1 (1+UX^2) +  2 UX UY P2
%   Pi_yy = -2 P1 (1+UY^2) +  2 UX UY P2
%   Pi_xy =  P2 ( 2 + UX^2 + UY^2 )

    A = params.A;
    Nx = params.Nx;
    Ny = params.Ny;
    Nt = Nx.*Ny;


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

    m11 = (-3).*T.^2.*(1+ux.^2+uy.^2).^(1/2);
    m12 = (1/2).*(1+ux.^2+uy.^2).^(-1/2).*((-3).*T.^3.*ux+(-4).*P1.*ux.*(1+uy.^2)+(-2).*P2.*uy.*(2+(-1).*ux.^2+uy.^2));
    m13 = (1/2).*(1+ux.^2+uy.^2).^(-1/2).*(((-3).*T.^3+4.*P1.*(1+ux.^2)).*uy+(-2).*P2.*ux.*(2+ux.^2+(-1).*uy.^2));

    m11 = (-3).*T.^2.*(1+ux.^2+uy.^2).^(1/2);
    m12 = (-3/2).*T.^3.*ux.*(1+ux.^2+uy.^2).^(-1/2);
    m13 = (-3/2).*T.^3.*uy.*(1+ux.^2+uy.^2).^(-1/2);

    m22 = (3/4).*T.^3.*(1+ux.^2+uy.^2).^(-1/2).*(2+ux.^2+2.*uy.^2);
    m23 = (-3/4).*T.^3.*ux.*uy.*(1+ux.^2+uy.^2).^(-1/2);

    m32 = (-3/4).*T.^3.*ux.*uy.*(1+ux.^2+uy.^2).^(-1/2);
    m33 = (3/4).*T.^3.*(1+ux.^2+uy.^2).^(-1/2).*(2+2.*ux.^2+uy.^2);


    R1 = (-3).*T.^2.*Tdx.*ux+(-3/2).*T.^3.*uxdx+(-3).*T.^2.*Tdy.*uy+(-3/2).*T.^3.*uydy;
    R2 = (3/2).*T.^2.*Tdx+(3/4).*T.^3.*ux.*uxdx+(3/2).*T.^3.*uxdy.*uy+(-3/4).*T.^3.*ux.*uydy;
    R3 = (3/2).*T.^2.*Tdy+(-3/4).*T.^3.*uxdx.*uy+(3/2).*T.^3.*ux.*uydx+(3/4).*T.^3.*uy.*uydy;
    R4 = P1.*(2+ux.^2+uy.^2);
    R5 = P2.*(2+ux.^2+uy.^2);

    mdet = (-1).*m11.*m23.*m32+m11.*m22.*m33;

    i11 = (-1).*m23.*m32+m22.*m33;
    i12 = m13.*m32+(-1).*m12.*m33;
    i13 = (-1).*m13.*m22+m12.*m23;

    i22 = m11.*m33;
    i23 = (-1).*m11.*m23;

    i32 = (-1).*m11.*m32;
    i33 = m11.*m22;



    i11 = i11./mdet;
    i12 = i12./mdet;
    i13 = i13./mdet;

    i22 = i22./mdet;
    i23 = i23./mdet;

    i32 = i32./mdet;
    i33 = i33./mdet;

    Tdt  = -(i11.*R1+i12.*R2+i13.*R3);
    uxdt = -(i22.*R2+i23.*R3);
    uydt = -(i32.*R2+i33.*R3);
    P1dt = 0.*i11;
    P2dt = 0.*i11;


    udt = zeros(Nx.*Ny.*5, 1);

    udt(0*Nt+1 : 1*Nt) = Tdt(:);
    udt(1*Nt+1 : 2*Nt) = uxdt(:);
    udt(2*Nt+1 : 3*Nt) = uydt(:);
    udt(3*Nt+1 : 4*Nt) = P1dt(:);
    udt(4*Nt+1 : 5*Nt) = P2dt(:);


end
