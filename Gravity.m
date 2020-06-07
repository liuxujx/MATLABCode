function [outp]  =  Gravity(inp,const)
% �����������ٶ�
% ����״̬��������
% ����������ٶ���λ��ʸ����γ�ȷ����ϵķ���
gE = const.gE; %��������������ٶ� [m/s^2]
GM = const.GM; %�������� [m^3/s^2]  
R0 = const.R0; %����������ƽ���뾶 [m]
g0 = const.g0; %�����������ٶ� [m/s^2]
a0 = const.a0; %���ƽ���뾶 [m]
J2 = const.J2; %2�״�гϵ�� �޵�λ

Rscale = const.Rscale;
Vscale = const.Vscale;
ascale = const.ascale;
tscale = const.tscale;

scale = const.scale; %�Ƿ������ٻ���־
if scale == 0
    r = inp(1);
    theta = inp(2);
    phi = inp(3);
    V = inp(4);
    gamma = inp(5);
    psi = inp(6);
elseif scale == 1 %��ԭΪ�����ٱ���
    r = inp(1)*Rscale;
    theta = inp(2);
    phi = inp(3);
    V = inp(4)*Vscale;
    gamma = inp(5);
    psi = inp(6);
else
    error('�����ٻ���־ scale=0 / 1');
end

%��������������λ��ʸ�������ϵĴ�С��û�з���
gr_i = GM/r^2;
%J2��������λ��ʸ�������ϵĴ�С��û�з���
gr_J2 = GM/r^2*(3/2*J2*(a0/r)^2*(1-3*sin(phi)^2));
%J2��������γ�ȷ����ϵĴ�С��û�з���
gphi_J2 = GM/r^2*(3*J2*(a0/r)^2*sin(phi)*cos(phi));

gr = gr_i+gr_J2;
gphi = gphi_J2;

outp = {gr,gphi};
end

