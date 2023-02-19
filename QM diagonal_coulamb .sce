//07/09/2022
clc
clear
n=1001//input("enter number of datapoints : ")
h=6.63D-34
Factor = 1//input("Input the Mass Factor : ")
//m=207*9.1D-31 (for muon, m = 207*me)
//e=1.6D-19
//for helium Z = 2 r_min = 1D-15 r_max = 1D-9
rmin=1D-6//input("enter r_min : ")
//xmin = input("enter the value of xmin : ")
//xmax = input("enter the value of xmax : ")
rmax=1//input("enter r_max : ")
ze = input("Enter the atomic number : ")
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
V=zeros(1,n)
V=-((E0*e)./r)/1D-9
A=eye(n-2,n-2)
v=diag(V(2:n-1))
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
//energy from Bohr spectrum
for i=1:counter
    Bohr(i)=-13.6/(i*i)
end
disp("Bound state eigen values from Bohr spectrum : ",Bohr)
for i=1:counter//percentage error
 err(i)=(abs(Bohr(i)-Z(i))/abs(Bohr(i)))*100
end
disp("Percentage error in eigen values ",err)
//--------------------------Normalizations & Expectation values-------------------------
vnew = V(2:n-1)
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
figure
for i=1:counter
    subplot(2,2,i)
    plot(rnew*1D9,u_r(:,i))
    xgrid()
    //m=gca()
    //m.x_ticks=tlist(["ticks","locations","labels"],[0.000D-10,5.000D-10,1.000D-9],["0e00","5e-10","1e-9"])
    title("U(r) vs r plot ","color","black","Fontsize","5","Fontname","2")
    xlabel("r----->","color","brown","Fontsize","3")
    ylabel("U(r)----->","color","brown","Fontsize","3")
end
figure
for i=1:counter
    subplot(2,2,i)
    plot(rnew*1D9,neigen(:,i))
    //m=gca()
    //m.x_ticks=tlist(["ticks","locations","labels"],[0.500D-10,5.000D-10,1.000D-9],["0.500e-10","5e-10","1e-9"])
    title("|U(r)|^2 plotting for bound state"+string(i),"color","black","Fontsize","5","Fontname","2")
    xlabel("r----->","color","brown","Fontsize","3")
    ylabel("|U(r)|^2----->","color","brown","Fontsize","3")
    xgrid()
end
//----------Radial R(r)------------------
for i=1:counter
    R(:,i)=u_r(:,i)./rnew
end
figure
for i=1:counter
    subplot(2,2,i)
    plot(rnew*1D9,R(:,i))
    xgrid()
    //m=gca()
    //m.x_ticks=tlist(["ticks","locations","labels"],[0.000D-10,5.000D-10,1.000D-9],["0e00","5e-10","1e-9"])
    title('R(r) vs r Plot','Fontsize','5','Fontname',2)
    xlabel('r-------->>','color','brown','Fontsize','4','Fontname',2)
    ylabel('R(r)= u(r)/r ---->>','color','brown','Fontsize','4','Fontname',3)

end
//------------Compare Radial Part--------------------------
a=0.529*1D-10
x=a^(-3/2)
r_a=rnew./a
R10 = 2*x*exp(-r_a)
R20 = (1/sqrt(2))*x*(1-((1/2).*r_a)).*exp(-r_a/2)
R30 = (2/sqrt(27))*x*(1-((2/3).*r_a)+((2/27).*((r_a)^2))).*exp(-r_a/3)
show_window(4)
for i=1:counter
    subplot(1,3,1)
    plot(rnew,R10,'-*r')
    plot(rnew,R(:,1),'-')
    legend("R_theo","R(r)")
    m=gca()
    m.x_ticks=tlist(["ticks","locations","labels"],[0.000D-10,5.000D-10,1.000D-9],["0e00","5e-10","1e-9"])
    title('Comparision of R_theo & R(r)','Fontsize','5','Fontname',2)
    xlabel('r----->','color','brown','Fontsize','4','Fontname',3)
    ylabel('R_theo & R(r)---->','color','brown','Fontsize','4','Fontname',3)
    xgrid()
    subplot(1,3,2)
    plot(rnew,R20,'-*r')
    plot(rnew,-R(:,2),'-')
    legend("R_theo","R(r)")
    m=gca()
    m.x_ticks=tlist(["ticks","locations","labels"],[0.000D-10,5.000D-10,1.000D-9],["0e00","5e-10","1e-9"])
    title('Comparision of R_theo & R(r)','Fontsize','5','Fontname',2)
    xlabel('r----->','color','brown','Fontsize','4','Fontname',3)
    ylabel('R_theo & R(r)---->','color','brown','Fontsize','4','Fontname',3)
    xgrid()
    subplot(1,3,3)
    plot(rnew,R30,'-*r')
    plot(rnew,R(:,3),'-')
    legend("R_theo","R(r)")
    m=gca()
    m.x_ticks=tlist(["ticks","locations","labels"],[0.000D-10,5.000D-10,1.000D-9],["0e00","5e-10","1e-9"])
    title('Comparision of R_theo & R(r)','Fontsize','5','Fontname',2)
    xlabel('r----->','color','brown','Fontsize','4','Fontname',3)
    ylabel('R_theo & R(r)---->','color','brown','Fontsize','4','Fontname',3)
    xgrid()
end
//-----------R_MP & Orthogonality--------------------------
for i=1:counter
    [MP IND]=max(neigen(:,i))
    r_MP(i) = rnew(IND)
end
disp("The Most probable value of r_MP at which u_max: ",r_MP)
R=0.528*1D-10//Bohr radius
for j=1:counter 
    B_R(j)=R*j*j
end
disp("The value of R_MP from theoretical formula : ",B_R)
//Orthogonality checking-------------------------------------
for i = 1:counter
    for j=1:counter
        if i~=j then
            neigen4 = u_r(:,i).*u_r(:,j)
            ortho(i,j)=inttrap(rnew,neigen4)
        end
    end
end
disp("Orthogonality check : ",ortho)




