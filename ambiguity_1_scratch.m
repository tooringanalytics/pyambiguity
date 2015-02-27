

function ambiguity_1_scratch (u_basic, fcode, f_basic, F, K, T, N, sr, plotdir, signal_name)

    mkdir(plotdir);

    %dmpdat('u_basic', u_basic);

    %if fcode == 1
    %    dmpdat('f_basic', f_basic);
    %end
    %hbrk("dmp");
    savdat(signal_name, 'u_basic', u_basic);
    m_basic=length(u_basic); %returns the largest array dimension in matrix
    df=F/K/m_basic;
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
        savdat(signal_name, 'ud', ud);

        ao=ones(r,m_basic);
        savdat(signal_name, 'ao', ao);

        m=m_basic*r;

        savdat(signal_name, 'ao_dot_ud', ao*ud);

        u_basic=reshape(ao*ud,1,m);
        savdat(signal_name, 'u_basic_reshaped', u_basic);

        uamp=abs(u_basic);
        savdat(signal_name, 'uamp', uamp);

        phas=angle(u_basic);
        savdat(signal_name, 'phas', phas)

        u=u_basic;
        if fcode==1
            ff=diag(f_basic);
            savdat(signal_name, 'ff', ff);
            vecprod = ao*ff;
            savdat(signal_name, 'vecprod', vecprod);
            vecprod_reshaped = reshape(ao*ff, 1, m);
            savdat(signal_name, 'vecprod_reshaped', vecprod_reshaped);
            cumsummed = cumsum(vecprod_reshaped);
            savdat(signal_name, 'cumsummed', cumsummed);
            add_term = 2*pi*dt*cumsummed;
            savdat(signal_name, 'add_term', add_term);
            phas=2*pi*dt*cumsum(reshape(ao*ff,1,m))+phas;
            savdat(signal_name, 'phas_add', phas);
            comprod=j*phas;
            savdat(signal_name, 'comprod', comprod);
            uexp=exp(j*phas);
            savdat(signal_name, 'uexp', uexp);
            u=uamp.*uexp;
            savdat(signal_name, 'u', u);
        end;
    end

    t=[0:r*m_basic-1]/r;
    savdat(signal_name, 't', t);

    tscale1=[0 0:r*m_basic-1 r*m_basic-1]/r;
    savdat(signal_name, 'tscale1', tscale1);

    dphas=[NaN diff(phas)]*r/2/pi;
    savdat(signal_name, 'dphas', dphas);

    %brk('[0 abs(uamp) 0]', [0 abs(uamp) 0]);
    figure(1), clf, hold off;
    subplot(3,1,1);
    savdat(signal_name, 'fig1plot311', [0 abs(uamp) 0]);
    plot(tscale1,[0 abs(uamp) 0],'linewidth',1.5);
    ylabel(' Amplitude ');
    axis([-Inf Inf 0 1.2*max(abs(uamp))]);

    subplot(3,1,2);
    savdat(signal_name, 'fig1plot312', phas);
    plot(t, phas,'linewidth',1.5);
    axis([-Inf Inf -Inf Inf]);
    ylabel(' Phase [rad] ');

    subplot(3,1,3);
    savdat(signal_name, 'fig1plot312', dphas*ceil(max(t)));
    plot(t,dphas*ceil(max(t)),'linewidth',1.5);
    axis([-Inf Inf -Inf Inf]);
    xlabel(' \itt / t_b ');
    ylabel(' \itf * Mt_b ');

    print('-dsvg', [ plotdir '/' signal_name '_fig_1' '.svg']);
    dtau=ceil(T*m)*dt/N;
    savdat(signal_name, 'dtau', dtau);
    tau=round([0:1:N]*dtau/dt)*dt;
    savdat(signal_name, 'tau', tau);


    f=[0:1:K]*df;
    f=[-fliplr(f) f];
    savdat(signal_name, 'f', f);
    %brk('Tm', ceil(T*m))
    %brk('m_plus_Tm', m+ceil(T*m))
    %brk('u.T', u');
    mat1=spdiags(u',0,m+ceil(T*m),m);
    savdat(signal_name, 'mat1', mat1);
    %brk('mat1', mat1)
    u_padded=[zeros(1,ceil(T*m)),u,zeros(1,ceil(T*m))];
    savdat(signal_name, 'u_padded', u_padded);
    cidx=[1:m+ceil(T*m)];
    savdat(signal_name, 'cidx', cidx);
    %brk('cidx', cidx)
    ridx=round(tau/dt)';
    savdat(signal_name, 'ridx', ridx);
    %brk('ar1', cidx(ones(N+1,1),:));
    %brk('ar2', ridx(:,ones(1,m+ceil(T*m))));
    index = cidx(ones(N+1,1),:) + ridx(:,ones(1,m+ceil(T*m)));
    savdat(signal_name, 'index', index);
    %brk('index', index);
    %brk('u_padded', u_padded)
    %brk('u_padded(index)', u_padded(index))
    mat2 = sparse(u_padded(index));
    savdat(signal_name, 'mat2', mat2);
    %brk('mat2', mat2)
    uu_pos=mat2*mat1;
    savdat(signal_name, 'uu_pos', uu_pos);
    %brk('uu_pos', uu_pos);
    %clear mat2 mat1
    e=exp(-j*2*pi*f'*t);
    savdat(signal_name, 'e', e);
    %brk('e', e);
    %brk('uu_pos', uu_pos)
    %brk('uu_pos_dash', uu_pos');
    %brk('e*uu_pos.conj()', e*uu_pos')
    a_pos=abs(e*uu_pos');
    savdat(signal_name, 'a_pos', a_pos);
    %brk('a_pos(121, 60)', a_pos(121, 60));
    %brk('a_pos', a_pos);
    a_pos=a_pos/max(max(a_pos));
    savdat(signal_name, 'a_pos_norm', a_pos);

    %brk('a_pos', a_pos);
    a=[flipud(conj(a_pos(1:K+1,:))) fliplr(a_pos(K+2:2*K+2,:))];
    savdat(signal_name, 'a', a);
    %brk('a(1, 113)', a(1, 113));
    delay=[-fliplr(tau) tau];
    savdat(signal_name, 'delay', delay);
    %brk('delay', delay);
    freq=f(K+2:2*K+2)*ceil(max(t));
    savdat(signal_name, 'freq', freq);
    %brk('f(K+2:2*K+2)', f(K+2:2*K+2))
    %brk('maxT', ceil(max(t)))
    %brk('freq', freq);
    delay=[delay(1:N) delay(N+2:2*N)];
    savdat(signal_name, 'delay_repl', delay);
    %brk('delay', delay);
    a=a(:,[1:N,N+2:2*N]);
    savdat(signal_name, 'a_repl', a);
    %brk('a', a);
    [amf amt]=size(a);
    savdat(signal_name, 'amf_amt', [amf amt]);
    %brk('amf', amf);
    %brk('amt', amt);
    cm=zeros(64,3);
    cm(:,3)=ones(64,1);
    %brk('cm', cm);

    figure(2), clf, hold off;
    %brk('[zeros(1,amt);a]', [zeros(1,amt);a]);
    savdat(signal_name, 'fig2meshx', delay);
    savdat(signal_name, 'fig2meshy', [0 freq]);
    savdat(signal_name, 'fig2meshz', [zeros(1,amt);a]);
    mesh(delay, [0 freq], [zeros(1,amt);a]);
    hold on

    savdat(signal_name, 'fig2surfx', delay);
    savdat(signal_name, 'fig2surfy', [0 0]);
    savdat(signal_name, 'fig2surfz', [zeros(1,amt);a(1,:)]);
    surface(delay, [0 0], [zeros(1,amt);a(1,:)]);
    colormap(cm);

    view(-40,50);
    axis([-Inf Inf -Inf Inf 0 1]);
    xlabel(' {\it\tau}/{\itt_b}','Fontsize',12);
    ylabel(' {\it\nu}*{\itMt_b}','Fontsize',12);
    zlabel(' |{\it\\chi}({\it\tau},{\it\nu})| ','Fontsize',12);
    hold off;
    print('-dsvg', [ plotdir '/' signal_name '_fig_2' '.svg']);

    savdat(signal_name, 'delay_final', delay);

    savdat(signal_name, 'freq_final', freq);

    savdat(signal_name, 'a_final', a);

endfunction

function dmpdat (nm, ex)
    printf('%s\n', nm);
    disp(ex);
    printf('size(%s)\n', nm);
    disp(size(ex));
end

function savdat(signal_name, exname, ex)

    odir = ['check_data' '/' signal_name];
    mkdir(odir);
    save('-mat', [odir '/' exname '.' 'mat'], 'ex')
    %hbrk('');

end

function hbrk (msg)

    printf('%s\n', msg);
    quit();

end

function brk (nm, ex)

    dmpdat(nm, ex);
    quit();

end

