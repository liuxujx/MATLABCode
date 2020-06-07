function [V_sound]= Vsound(altitude,planetary)
%����������٣�����߶ȣ��������
% ���򣺻��� 1976 US ��׼����ģ���������
% ���ǣ��ο����ף�Advances in sapcecraft atmospheric entry gudiance, P50-51
% ����=altitude �߶� m (0=< alt =< 85000 m)
% ���=Vsound ���� m/s
%
if planetary==1  % 1=���� 2=����
    alt = altitude/1000; %[km]
    if(alt > 85.0)
        alt = 85.0;
    end
    if(alt < 0.0)
        alt = 0.0;
    end
    a = 340.29112;
    b = -0.87245462;
    c = -300.10255;
    d = 0.3077688;
    e = 106.65518;
    f = -0.052091446;
    g = -18.273173;
    c_h = 0.004143208;
    c_i = 1.485658;
    c_j = -0.00012124848;
    c_k = -0.045201417;
    h_s = sqrt(alt);
    dum1 = a + c*h_s + e*alt + g*h_s*h_s*h_s + c_i*alt*alt + c_k*power(h_s,5.0);
    dum2 = 1.0+b*h_s + d*alt + f*h_s*h_s*h_s + c_h*alt*alt + c_j*power(h_s,5.0);
    V_sound=dum1/dum2;
elseif planetary==2
    alt = altitude; %[m]
    if alt<=10e3 || alt>=120e3
        warning('���Ǵ�������ģ����[10~120]km֮�ⲻ׼ȷ');
    end
    b0 = 223.8;
    b1 = -0.0002004;
    b2 = -1.5888e-8;
    b3 = 1.404e-13;
    V_sound = b0+b1*alt+b2*alt.^2+b3*alt.^3;
end