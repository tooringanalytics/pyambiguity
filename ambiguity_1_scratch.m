
function brk (nm, ex)
    printf('%s\n', nm);
    disp(ex);
    printf('size(%s)\n', nm);
    disp(size(ex));
    quit();
endfunction

% clear all
% t=0:1:255;
% st=chirp(t,0,10,100);
% wt=kaiser(256,5)';
% st1=st.*wt;
st1=ones(1,51);
%st1=1;
u_basic= st1; %input( 'Signal elements (row complex vector, each element last tb sec) = ? ');
m_basic=length(u_basic); %returns the largest array dimension in matrix
%fcode=input(' Allow frequency coding (yes=1, no=0) = ? ');
%if fcode==1
%    f_basic=input(' Frequency coding in units of 1/tb (row vector of same length) = ? ');
%end

fcode=1
f_basic=0.0031 * (-25:25)
F=6
K=60
T=1.1
N=60
sr=10

%F=input(' Maximal Doppler shift for ambiguity plot [in units of 1/Mtb] (e.g., 1)= ? ');
%K=input(' Number of Doppler grid points for calculation (e.g., 100) = ? ');
%F=float(F);
df=F/K/m_basic;
%T=input(' Maximal Delay for ambiguity plot [in units of Mtb] (e.g., 1)= ? ');
%N=input(' Number of delay grid points on each side (e.g. 100) = ? ');
%sr=input(' Over sampling ratio (>=1) (e.g. 10)= ? ');
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
        comprod=j*phas;
        uexp=exp(j*phas);
        u=uamp.*uexp;
        %brk("u", u);
    end
end

t=[0:r*m_basic-1]/r;
tscale1=[0 0:r*m_basic-1 r*m_basic-1]/r;
dphas=[NaN diff(phas)]*r/2/pi;

figure(1), clf, hold off;
subplot(3,1,1);
plot(tscale1,[0 abs(uamp) 0],'linewidth',1.5);
ylabel(' Amplitude ');
axis([-inf inf 0 1.2*max(abs(uamp))]);

subplot(3,1,2);
plot(t, phas,'linewidth',1.5);
axis([-inf inf -inf inf]);
ylabel(' Phase [rad] ');

subplot(3,1,3);
plot(t,dphas*ceil(max(t)),'linewidth',1.5);
axis([-inf inf -inf inf]);
xlabel(' \itt / t_b ');
ylabel(' \itf * Mt_b ');


dtau=ceil(T*m)*dt/N;
tau=round([0:1:N]*dtau/dt)*dt;


f=[0:1:K]*df;
f=[-fliplr(f) f];
%brk('Tm', ceil(T*m))
%brk('m_plus_Tm', m+ceil(T*m))
%brk('u.T', u');
mat1=spdiags(u',0,m+ceil(T*m),m);
%brk('mat1', mat1)
u_padded=[zeros(1,ceil(T*m)),u,zeros(1,ceil(T*m))];
cidx=[1:m+ceil(T*m)];
%brk('cidx', cidx)
ridx=round(tau/dt)';
%brk('ar1', cidx(ones(N+1,1),:));
%brk('ar2', ridx(:,ones(1,m+ceil(T*m))));
index = cidx(ones(N+1,1),:) + ridx(:,ones(1,m+ceil(T*m)));
%brk('index', index);
%brk('u_padded', u_padded)
%brk('u_padded(index)', u_padded(index))
mat2 = sparse(u_padded(index));
%brk('mat2', mat2)
uu_pos=mat2*mat1;
%brk('uu_pos', uu_pos);
%clear mat2 mat1
e=exp(-j*2*pi*f'*t);
%brk('e', e);
%brk('uu_pos', uu_pos)
%brk('uu_pos_dash', uu_pos');
%brk('e*uu_pos.conj()', e*uu_pos')
a_pos=abs(e*uu_pos');
%brk('a_pos(121, 60)', a_pos(121, 60));
%brk('a_pos', a_pos);
a_pos=a_pos/max(max(a_pos));
%brk('a_pos', a_pos);
a=[flipud(conj(a_pos(1:K+1,:))) fliplr(a_pos(K+2:2*K+2,:))];
%brk('a(1, 113)', a(1, 113));
delay=[-fliplr(tau) tau];
%brk('delay', delay);
freq=f(K+2:2*K+2)*ceil(max(t));
%brk('f(K+2:2*K+2)', f(K+2:2*K+2))
%brk('maxT', ceil(max(t)))
%brk('freq', freq);
delay=[delay(1:N) delay(N+2:2*N)];
%brk('delay', delay);
a=a(:,[1:N,N+2:2*N]);
%brk('a', a);
[amf amt]=size(a);
%brk('amf', amf);
%brk('amt', amt);
cm=zeros(64,3);
cm(:,3)=ones(64,1);
%brk('cm', cm)

figure(2), clf, hold off
%brk('[zeros(1,amt);a]', [zeros(1,amt);a])
mesh(delay, [0 freq], [zeros(1,amt);a])
hold on

surface(delay, [0 0], [zeros(1,amt);a(1,:)])
colormap(cm)

view(-40,50)
axis([-inf inf -inf inf 0 1])
xlabel(' {\\it\\tau}/{\\itt_b}','Fontsize',12);
ylabel(' {\\it\\nu}*{\\itMt_b}','Fontsize',12);
zlabel(' |{\\it\\chi}({\\it\\tau},{\\it\\nu})| ','Fontsize',12);
hold off

