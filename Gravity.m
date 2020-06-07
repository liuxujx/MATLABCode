function [gr, gphi] = Gravity(state,auxdata)
% �����������ٶ�
% ����״̬��������
% ����������ٶ���λ��ʸ����γ�ȷ����ϵĴ�С������ο����� P43 ͼ2-5
% �ο����ף�Զ�̻�����ж���ѧ���Ƶ���P40-48
gE = auxdata.gE; %��������������ٶ� [m/s^2]
GM = auxdata.GM; %�������� [m^3/s^2]  
R0 = auxdata.R0; %����������ƽ���뾶 [m]
g0 = auxdata.g0; %�����������ٶ� [m/s^2]
a0 = auxdata.a0; %���ƽ���뾶 [m]
J2 = auxdata.J2; %2�״�гϵ�� �޵�λ

scale = auxdata.scale; %�Ƿ������ٻ���־  0=������ 1=������
Rscale = auxdata.Rscale;
Vscale = auxdata.Vscale;
ascale = auxdata.ascale;
tscale = auxdata.tscale;

r = state(1);
% theta = state(2);
phi = state(3);
% V = state(4);
% gamma = state(5);
% psi = state(6);

if scale == 0 %�������������ٶ�
    %��������������λ��ʸ�������ϵĴ�С��û�з���
    gr_i = GM/r^2;
    %J2��������λ��ʸ�������ϵĴ�С��û�з���
    gr_J2 = GM/r^2*(3/2*J2*(a0/r)^2*(1-3*sin(phi)^2));
    %J2��������γ�ȷ����ϵĴ�С��û�з���
    gphi_J2 = GM/r^2*(3*J2*(a0/r)^2*sin(phi)*cos(phi));
elseif scale == 1 %�����ٻ��������ٶ�
    %��������������λ��ʸ�������ϵĴ�С��û�з���
    gr_i = 1/r^2;
    %J2��������λ��ʸ�������ϵĴ�С��û�з���
    gr_J2 = 1/r^2*(3/2*J2*(a0/Rscale/r)^2*(1-3*sin(phi)^2));
    %J2��������γ�ȷ����ϵĴ�С��û�з���
    gphi_J2 = 1/r^2*(3*J2*(a0/Rscale/r)^2*sin(phi)*cos(phi));
end

gr = gr_i+gr_J2;
gphi = gphi_J2;

end

