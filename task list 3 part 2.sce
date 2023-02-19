clc
clear
n=input("enter number of datapoints : ")
h=6.63D-34
Factor = 1
//m=207*9.1D-31 (for muon, m = 207*me)
//e=1.6D-19
//for helium Z = 2 r_min = 1D-15 r_max = 1D-9
rmin=1D-6//input("enter r_min : ")
//xmin = input("enter the value of xmin : ")
//xmax = input("enter the value of xmax : ")
rmax=1//input("enter r_max : ")
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
a0= input("Enter the value of a0 : ")
V =-(((E0*e)*exp(-new/a0))./new)
A=eye(n-2,n-2)
v=diag(V)
//disp(A)
//----------------------By inbuilt command-------------------------
//tic()
D=(-2*C)*ones(n-2,1)
A1=diag(D)
//disp(A1)
C1=C*ones(n-3,1)
A2=diag(C1,1)
//disp(A2)
A3=diag(C1,-1)
//disp(A3)
//Hermitian matrix (tridiagonal)
H=A1+A2+A3+v
//disp("By inbuilt : ",H)
[a1 a2]=spec(H)
Z=spec(H)
//disp("Corresponding eigen vectors : ",a1)
//disp("Eigen values : ",a2)
counter = 0
for i=1:n-2
    if a2(i,i)<0
        counter=counter+1
    end
end
disp("No. of Bound states : "+string(counter))
disp("Eigen values for Bound state : ",Z(1:counter))
//--------------------------Normalizations & Expectation values-------------------------
vnew = V
for i = 1:counter
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(rnew,eigenvec)
    u_r(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = u_r(:,i).*u_r(:,i)
    neigen1(:,i) = rnew.*neigen(:,i)
    neigen2(:,i) = rnew.*neigen1(:,i)
    neigen3(:,i) = vnew'.*neigen(:,i)
    ex_r(i)=inttrap(rnew,neigen1(:,i))
    ex_r2(i)= inttrap(rnew,neigen2(:,i))
    ex_V(i) = inttrap(rnew,neigen3(:,i))
end
for i = 1:n-3
    for j = 1:counter
        r2(i) = (rnew(i)+rnew(i+1))/2
        mid_u(i,j) = (u_r(i,j)+u_r(i+1,j))/2
        diff_u(i,j) = (u_r(i+1,j)-u_r(i,j))/dr
    end
end

for i = 1:n-4
    for j = 1:counter
        r3(i) = (r2(i) + r2(i+1))/2
        mid2_u(i,j) = (u_r(i+2,j) + 2*u_r(i+1,j)+ u_r(i,j))/4
        diff2_u(i,j) = (u_r(i+2,j) - 2*u_r(i+1,j) + u_r(i,j))/(dr*dr)
    end
end
//----------------Uncertainty check--------------------------------------
Un = (4*%pi)/h

for i=1:counter
    y2(:,i) = mid_u(:,i).*diff_u(:,i)
    y3(:,i) = mid2_u(:,i).*diff2_u(:,i)
    ex_p(i) = -1*%i*(h/(2*%pi))*inttrap(r2,y2(:,i))
    ex_p2(i) = -1*(h/(2*%pi))**2*inttrap(r3,y3(:,i))
    sig_r(i) = sqrt(ex_r2(i) - (ex_r(i)*ex_r(i)))
    sig_p(i) = sqrt(ex_p2(i) - (ex_p(i)*ex_p(i)))
    ex_K(i) = ex_p2(i)/(2*m*(e/ze))
    un(i) = Un*(sig_r(i).* sig_p(i))
end
//-----------------total energy------------------------
E = ex_V + ex_K
disp("Expectation value <r> =",ex_r)
disp("Expectation value <r2> =",ex_r2)
disp("Expectation value <p> =",ex_p)
disp("Expectation value <p2> =",ex_p2)
disp("Expectation value <V> =",ex_V)
disp("Expectation value <KE> =",ex_K)
disp("Total Energy for all bound state values <V>+<KE> = ",E)
disp("Standard deviation of r =",sig_r)
disp("Standard deviation of p =",sig_p)
disp("Uncertanity Product (hbar/2) = ",un)
for i=1:counter
    if abs(E(i)-Z(i))<=0.1
    disp("Energy Eigenvalues Correct for boundstate"+string(i))
    else
    disp("Energy Eigenvalues Incorrect for boundstate"+string(i))
    end
end
show_window(2)
for i=1:counter
    subplot(1,counter,i)
    plot(rnew*1D9,u_r(:,i),'m','linewidth',2)
    title("<<u_r vs r(in nm) Plot for bound state "+string(i)+">>",'color','brown','font',2,'fontsize',2)
    xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
    ylabel("u_r----------->",'color','brown','font',2,'fontsize',4)
    xgrid()
end
//----------Radial R(r)------------------
for i=1:counter
    R(:,i)=u_r(:,i)./rnew
end
show_window(3)
for i=1:counter
    subplot(1,counter,i)
    plot(rnew*1D9,neigen(:,i),'k','linewidth',2)
    title("<<|u_r|^2 vs r(in nm) Plot for bound state "+string(i)+">>",'color','brown','font',2,'fontsize',2)
    xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
    ylabel("|u_r|^2----------->",'color','brown','font',2,'fontsize',4)
    xgrid()
end
show_window(4)
for i=1:counter
    subplot(1,counter,i)
    plot(rnew*1D9,R(:,i),'r','linewidth',2)
    title("<<R_r vs r(in nm) Plot for bound state "+string(i)+">>",'color','brown','font',2,'fontsize',2)
    xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
    ylabel("R_r----------->",'color','brown','font',2,'fontsize',4)
    xgrid()
end


