clear all
set(0,'DefaultFigureVisible','on')

%Constants (all MKS, except energy which is in eV)
hbar=1.055e-34;m=9.110e-31;epsil=8.854e-12;q=1.602e-19;
%Variables
v0=10;
b=1e-9;
L0=2e-9;
%for scaling graph
if v0>1
    scalec=1;
else
    scalec=v0/v0*0.1;
end

%Lattice
Np=100;L=L0+2*b;a=L/(Np);

%graphical eigenvalues
En=(hbar^2)/(2*m*(a^2))/q;
%analytical eigenvalues
Ean=(((hbar*pi)^2)/(2*m*(L0^2))/q)*[1:Np].*[1:Np];

%hamiltonian matrix v0 determining potential difference outside the well
Z=zeros(1,Np);
for n=1:Np
    if n*a < b
        Z(n)=v0;
    elseif n*a > L-b
        Z(n)=v0;
    else
        Z(n)=0;
    end
end
Vm=diag(Z);

% Setting up hamiltonian matrix
T=(Vm*diag(ones(1,Np)))+(2*En*diag(ones(1,Np)))-(En*diag(ones(1,Np-1),1))-(En*diag(ones(1,Np-1),-1));
[V,D]=eig(T);D=diag(D);[Enum,ind]=sort(D);

f1=figure;
%plot wavefunction
hold on
%plotting the well
plot([-b 0], [v0 v0],'k');plot([L0 L0+b] ,[v0 v0],'k');plot ([0 L0] , [0 0 ],'k')
plot([0 0],[0 v0],'k');plot([L0 L0],[0 v0],'k');
x=linspace (-b,L0+b,100);
E=zeros(1,Np)';
Psi=zeros(1,Np);
plotpsi=zeros(100,10);
for i=1:Np
    E(i,1)=D(ind(i));
    if E(i,1)  > v0
        break
    end
    Psi=abs(V(:,ind(i)));
    plotpsi(:,i)=Psi;
%normalizing the finite well wavefunction and scaling
    Psi=Psi/norm(Psi);
    Psi=rescale(Psi,0,scalec);
    Psi=Psi-Psi(1,1)*ones(100,1);
    plot(x,Psi+E(i),'b');
end
% set(gca,'Fontsize',[25])
axis([-b (b+L0) 0 inf])
%xlabel('Eigenvalue Number , m --->');% Part (a)
%ylabel('E (eV) --->');% Part (a)
xlabel('Position');% Part (b)
ylabel(' Energy Level ');% Part (b)
f2=figure;
hold on
h1=plot(Enum,'bx');% Part (a)
h=plot(Ean,'b');% Part (a)
xlim([0 20]);
ylim([0 Enum(20)]);
xlabel('Eigenvalue Number , m --->');% Part (a)
ylabel('E (eV) --->');% Part (a)
grid on
