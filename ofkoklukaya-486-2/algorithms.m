clc
clear
close all
load hw2_data.mat;

%% Q1

%from IGRF outout

O1=zeros(3,25000);
O2=zeros(3,25000);
O3=zeros(3,25000);
B=zeros(3,25000);
A_oi=zeros(3,3,25000);
mag_orb=zeros(3,25000);
mdm_e=7.71*10^15; %The magnetic dipole moment of the Earth
inc= deg2rad(98.1245); %[rad] inclination
dist=7056198.5; %[m] distance
mu=3.986004418*10^14;
om_0=sqrt(mu/dist^3);
om_e=7.29*10^-5; %[rad/sec] the spin rate of the earth
time=1:25000;


for i=1:25000
    O3(:,i)=-pos_eci(:,i)/norm(pos_eci(:,i));
    O2(:,i)=(cross(-pos_eci(:,i),vel_eci(:,i)))/norm(cross(-pos_eci(:,i),vel_eci(:,i)));
    O1(:,i)=cross(O2(:,i),O3(:,i));
    A_oi(:,:,i)=[O1(:,i),O2(:,i),O3(:,i)];
    mag_orb(:,i)=A_oi(:,:,i)*mag_eci(:,i);
end

%from dipole moment
delay = 4360;
for step=1:25000
    t = step + delay;
    B(1,step)=(mdm_e/dist^3)*(cos(om_0*t)*(cos(eps)*sin(inc)-sin(eps)*cos(inc)*...
        cos(om_e*t))-sin(om_0*t)*sin(eps)*sin(om_e*t));
    B(2,step)=-(mdm_e/dist^3)*(cos(eps)*cos(inc)+sin(eps)*sin(inc)*cos(om_e*t));
    B(3,step)=(2*mdm_e/dist^3)*(sin(om_0*t)*(cos(eps)*sin(inc)-sin(eps)*cos(inc)*...
        cos(om_e*t))-2*sin(om_0*t)*sin(eps)*sin(om_e*t));
end


figure;
subplot(3,2,1);
plot(time,mag_orb(1,:));
title("x-direction of Magnetic field vector in orbit frame from IGRF");
xlabel("t(s)");
subplot(3,2,3);
plot(time,mag_orb(2,:));
title("y-direction of Magnetic field vector in orbit frame from IGRF");
xlabel("t(s)");
subplot(3,2,5);
plot(time,mag_orb(3,:));
title("z-direction of Magnetic field vector in orbit frame from IGRF");
xlabel("t(s)");


subplot(3,2,2);
plot(time,B(1,:));
title("x-direction of Magnetic field vector in orbit frame from dipole model");
xlabel("t(s)");
subplot(3,2,4);
plot(time,B(2,:));
title("y-direction of Magnetic field vector in orbit frame from dipole model");
xlabel("t(s)");
subplot(3,2,6);
plot(time,B(3,:));
title("z-direction of Magnetic field vector in orbit frame from dipole model");
xlabel("t(s)");


%%

%% Q2

for s=1:25000
    sun_orb(:,s)=A_oi(:,:,s)*sun_eci(:,s);
end
figure;
subplot(3,1,1);
plot(time,sun_orb(1,:));
title("x-direction of sun direction vector in orbit frame");
xlabel("t(s)");
subplot(3,1,2);
plot(time,sun_orb(2,:));
title("y-direction of sun direction vector in orbit frame");
xlabel("t(s)");
subplot(3,1,3);
plot(time,sun_orb(3,:));
title("z-direction of sun direction vector in orbit frame");
xlabel("t(s)");

%%

%% Q3

E=zeros(3,25000); %roll,pitch,yaw
E(2,:)=5*sin(time(:)/500);

figure;
subplot(3,1,1);
plot(time,E(1,:));
title("roll angle wrt orbit frame");
xlabel("t(s)");
subplot(3,1,2);
plot(time,E(2,:));
title("pitch angle wrt orbit frame");
xlabel("t(s)");
subplot(3,1,3);
plot(time,E(3,:));
title("yaw angle wrt orbit frame");
xlabel("t(s)");

%%

%% Q4

E_1=transpose(E);
A_os=zeros(3,3,25000);

for e=1:25000
    A_os(:,:,e)=angle2dcm(E_1(e,1),E_1(e,2),E_1(e,3));
    %for magnetometer
    mag_body(:,e)=A_os(:,:,e)*mag_orb(:,e) + randn(3,1)*300;
    %for sun sensor
    sun_body(:,e)=A_os(:,:,e)*sun_orb(:,e) + randn(3,1)*0.02;
end

figure;
subplot(3,2,1);
plot(time,mag_body(1,:));
title("x-direction of Magnetic field vector in body frame with noise");
xlabel("t(s)");
subplot(3,2,3);
plot(time,mag_body(2,:));
title("y-direction of Magnetic field vector in body frame with noise");
xlabel("t(s)");
subplot(3,2,5);
plot(time,mag_body(3,:));
title("z-direction of Magnetic field vector in body frame with noise");
xlabel("t(s)");


subplot(3,2,2);
plot(time,sun_body(1,:));
title("x-direction of sun direction vector in body frame with noise");
xlabel("t(s)");
subplot(3,2,4);
plot(time,sun_body(2,:));
title("y-direction of sun direction vector in body frame with noise");
xlabel("t(s)");
subplot(3,2,6);
plot(time,sun_body(3,:));
title("z-direction of sun direction vector in body frame with noise");
xlabel("t(s)");


%%

%% Q5

%%