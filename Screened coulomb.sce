clc
clear
n=1001;//input("enter number of datapoints : ")
h=6.63D-34
Factor = 1
//m=207*9.1D-31 (for muon, m = 207*me)
//e=1.6D-19
//for helium Z = 2 r_min = 1D-15 r_max = 1D-9
rmin=1D-6
//xmin = input("enter the value of xmin : ")
//xmax = input("enter the value of xmax : ")
rmax=1
ze = 1
e=ze*1.6D-19 
m=Factor*9.1D-31
//vt=xt_2
//v0=-100
E0 = 9*1D9
r=linspace(rmin,rmax,n)//range
new = r(2:n-1)*1D-9
rnew = new'
dr=((rmax-rmin)/(n-1))*1D-9
C=-((h/(2*%pi))^2)/(2*m*(dr^2)*(e/ze))
//--------for coulomb potential-------------------------------------
VC=zeros(1,n)
VC =-((E0*e)./r)/1D-9 
//--------------for screened coulomb potential---------------------
a0=[0.3*10^(-9) 0.5*10^(-9) 0.7*10^(-9)]
for k=1:3
V =-(((E0*e)*exp(-new/a0(k)))./new)
V_new(:,k)=V
end
plot(rnew,VC(2:n-1),'*y')
plot(rnew,V_new(:,1),'r')
plot(rnew,V_new(:,2),'b')
plot(rnew,V_new(:,3),'g')
title("<<<Coulomb Potential vs Screened Coulomb Potential>>>>",'color','brown','font',2,'fontsize',4)
legend("Coulomb Potential in eV","SC Potential for a=0.3 nm","SC Potential for a=0.5 nm","SC Potential for a=0.7 nm")
xlabel("r--------------->",'color','brown','font',2,'fontsize',4)
ylabel("Potential in eV------------>",'color','brown','font',2,'fontsize',4)
xgrid()
z=gca();
z.data_bounds=[0,-80;2*10^(-10),0]
