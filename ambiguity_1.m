% clear all
% t=0:1:255;
% st=chirp(t,0,10,100);
% wt=kaiser(256,5)';
% st1=st.*wt;
st1=ones(1,51);
%st1=1;
u_basic= st1; %input( 'Signal elements (row complex vector, each element last tb sec) = ? ');
m_basic=length(u_basic); %returns the largest array dimension in matrix
fcode=input(' Allow frequency coding (yes=1, no=0) = ? ');
if fcode==1
    f_basic=input(' Frequency coding in units of 1/tb (row vector of same length) = ? ');
end

F=input(' Maximal Doppler shift for ambiguity plot [in units of 1/Mtb] (e.g., 1)= ? ');
K=input(' Number of Doppler grid points for calculation (e.g., 100) = ? ');
%F=float(F);
df=F/K/m_basic;
T=input(' Maximal Delay for ambiguity plot [in units of Mtb] (e.g., 1)= ? ');
N=input(' Number of delay grid points on each side (e.g. 100) = ? ');
sr=input(' Over sampling ratio (>=1) (e.g. 10)= ? ');
r=ceil(sr*(N+1)/T/m_basic);

if r==1
    dt=1;
    m=m_basic;
    uamp=abs(u_basic);
    phas=uamp*0;
    phas=angle(u_basic);
    if fcode==1
        phas=phas+2*pi*cumsum(f_basic);
    end
    uexp=exp(j*phas);
    u=uamp.*uexp;
else
    dt=1/r;
    ud=diag(u_basic);
    ao=ones(r,m_basic);
    m=m_basic*r;
    u_basic=reshape(ao*ud,1,m);
    uamp=abs(u_basic);
    phas=angle(u_basic);
    u=u_basic;
    if fcode==1
        ff=diag(f_basic);
        phas=2*pi*dt*cumsum(reshape(ao*ff,1,m))+phas;
        uexp=exp(j*phas);
        u=uamp.*uexp;
    end
end

t=[0:r*m_basic-1]/r;
tscale1=[0 0:r*m_basic-1 r*m_basic-1]/r;

dphas=[NaN diff(phas)]*r/2/pi;

figure(1), clf, hold off
subplot(3,1,1)
plot(tscale1,[0 abs(uamp) 0],'linewidth',1.5)
ylabel(' Amplitude ')
axis([-inf inf 0 1.2*max(abs(uamp))])

subplot(3,1,2)
plot(t, phas,'linewidth',1.5)
axis([-inf inf -inf inf])
ylabel(' Phase [rad] ')

subplot(3,1,3)
plot(t,dphas*ceil(max(t)),'linewidth',1.5)
axis([-inf inf -inf inf])
xlabel(' \itt / t_b ')
ylabel(' \itf * Mt_b ')

dtau=ceil(T*m)*dt/N;
tau=round([0:1:N]*dtau/dt)*dt;
f=[0:1:K]*df;
f=[-fliplr(f) f];
mat1=spdiags(u',0,m+ceil(T*m),m);
u_padded=[zeros(1,ceil(T*m)),u,zeros(1,ceil(T*m))];
cidx=[1:m+ceil(T*m)];
ridx=round(tau/dt)';
index = cidx(ones(N+1,1),:) + ridx(:,ones(1,m+ceil(T*m)));
mat2 = sparse(u_padded(index));
uu_pos=mat2*mat1;
clear mat2 mat1
e=exp(-j*2*pi*f'*t);
a_pos=abs(e*uu_pos');
a_pos=a_pos/max(max(a_pos));
a=[flipud(conj(a_pos(1:K+1,:))) fliplr(a_pos(K+2:2*K+2,:))];
delay=[-fliplr(tau) tau];
freq=f(K+2:2*K+2)*ceil(max(t));
delay=[delay(1:N) delay(N+2:2*N)];
a=a(:,[1:N,N+2:2*N]);
[amf amt]=size(a);
cm=zeros(64,3);
cm(:,3)=ones(64,1);
figure(2), clf, hold off
mesh(delay, [0 freq], [zeros(1,amt);a])
hold on
surface(delay, [0 0], [zeros(1,amt);a(1,:)])
colormap(cm)
view(-40,50)
axis([-inf inf -inf inf 0 1])
xlabel(' {\it\tau}/{\itt_b}','Fontsize',12);
ylabel(' {\it\nu}*{\itMt_b}','Fontsize',12);
zlabel(' |{\it\chi}({\it\tau},{\it\nu})| ','Fontsize',12);
hold off

