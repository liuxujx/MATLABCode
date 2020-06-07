function [outp]=EntryDyna3DoF(t,state,control,auxdata)
% 3-DoF����ѧ����

gE = auxdata.gE; %��������������ٶ� [m/s^2]
GM=auxdata.GM;
R0=auxdata.R0; %����������ƽ���뾶 [m]
g0=auxdata.g0; %�����������ٶ� [m/s^2]
a0 = auxdata.a0; %���ƽ���뾶 [m]
J2 = auxdata.J2; %2�״�гϵ�� �޵�λ

scale = auxdata.scale; %�Ƿ������ٻ���־ 0=������ 1=������
Rscale=auxdata.Rscale;
Vscale=auxdata.Vscale;
ascale=auxdata.ascale;
tscale=auxdata.tscale;

d2r=auxdata.d2r;
r2d=auxdata.r2d;

CL=auxdata.CL;
CD=auxdata.CD;
Sr=auxdata.Sr;
mass=auxdata.mass;

rho0 = auxdata.rho0;
hs = auxdata.hs;

r=state(1);
theta=state(2);
phi=state(3);
V=state(4);
gamma=state(5);
psi=state(6);

sigma=control(1);

if scale==0
    h=r-R0;
    rho=rho0*exp(-h/hs);
    %�������������ٶ�
    D=rho*V^2*Sr*CD/(2*mass);
    L=rho*V^2*Sr*CL/(2*mass);
    %�������ٶ�
    [gr,gphi]=Gravity(state,auxdata);
    %��ת���ٶ�
    OM=auxdata.OM;
elseif scale==1
    h=(r-1)*Rscale;
    rho=rho0*exp(-h/hs);
    %�����ٵ��������������ٶ�
    D=rho*(Vscale*V)^2*Sr*CD/(2*mass*g0);
    L=rho*(Vscale*V)^2*Sr*CL/(2*mass*g0);
    %�������ٶ�
    [gr,gphi]=Gravity(state,auxdata);
    %��ת���ٶ�
    OM=auxdata.OM*tscale;
end

rdot=V*sin(gamma);
thetadot=V*cos(gamma)*sin(psi)/(r*cos(phi));
phidot=V*cos(gamma)*cos(psi)/r;
Vdot= -D - sin(gamma)*gr - cos(gamma)*cos(psi)*gphi + ...
    OM^2*r*cos(phi)*(sin(gamma)*cos(phi)-cos(gamma)*sin(phi)*cos(psi));
gammadot=1/V*( L*cos(sigma) + (V^2/r-gr)*cos(gamma) + 2*OM*V*...
    cos(phi)*sin(psi) + OM^2*r*cos(phi)*(cos(gamma)*cos(phi)+...
    sin(gamma)*cos(psi)*sin(phi)) + sin(gamma)*cos(psi)*gphi );
psidot=1/V*( L*sin(sigma)/cos(gamma) + V^2*cos(gamma)*sin(psi)*tan(phi)/r...
    - 2*OM*V*(tan(gamma)*cos(psi)*cos(phi)-sin(phi)) + ...
    OM^2*r*sin(psi)*sin(phi)*cos(phi)/cos(gamma) + sin(psi)/cos(gamma)*gphi );

outp=[rdot thetadot phidot Vdot gammadot psidot].';
end