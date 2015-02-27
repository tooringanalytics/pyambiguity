
function test_ambiguity_1_scratch (sigcodes=[])

    plotdir = 'matlab_plots';
    mkdir(plotdir);

    % Plot the signals in the sigcodes array

    for sigcode = sigcodes,

        f_basic = 0;

        if sigcode == 0 || sigcode == 1

            signal_name='Pulse';
            u_basic = 1;
            fcode = 0;
            F=4;
            K=60;
            T=1.1;
            N=60;
            sr=10;

        end;

        if sigcode == 0 || sigcode == 2

            signal_name='LFM';
            u_basic = ones(1, 51);
            fcode=1;
            f_basic=0.0031 * (-25:25);
            F=6;
            K=60;
            T=1.1;
            N=60;
            sr=10;

            disp(u_basic);

        end;


        if sigcode == 0 || sigcode == 3

            signal_name='Weighted LFM';
            u_basic = sqrt(chebwin(51, 50))';
            fcode = 1;
            f_basic = .0031 * (-25:25);
            F=6;
            K=60;
            T=1.1;
            N=60;
            sr=10;

        end;



        if sigcode == 0 || sigcode == 4

            signal_name='Costas 7';
            u_basic = ones(1, 7);
            fcode = 1;
            f_basic = [4 7 1 6 5 2 3]
            F = 12;
            K = 60;
            T = 1.1;
            N = 80;
            sr = 20;

        end;


        if sigcode == 0 || sigcode == 5

            signal_name='Barker 13';
            u_basic = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
            fcode = 0;
            F = 10;
            K = 60;
            T = 1.1;
            N = 60;
            sr = 10;

        end;


        if sigcode == 0 || sigcode == 6

            signal_name='Frank 16';
            u_basic = [1 1 1 1 1 j -1 -j 1 -1 1 -1 1 -j -1 j];
            fcode = 0;
            F = 10;
            K = 60;
            T = 1.1;
            N = 60;
            sr = 10;

        end;


        if sigcode == 0 || sigcode == 7

            signal_name='P4 25';
            u_basic = exp(j*pi*(1/25*(0:24).^2-(0:24)));
            fcode = 0;
            F = 15;
            K = 80;
            T = 1.1;
            N = 80;
            sr = 20;

        end;


        if sigcode == 0 || sigcode == 8

            signal_name='Complementary Pair';
            u_basic = [1 1 -1 0 0 0 0 0 0 0 1 j 1];
            fcode = 0;
            F = 10;
            K = 60;
            T = 1.1;
            N = 60;
            sr = 10;

        end;


        if sigcode == 0 || sigcode == 9

            signal_name='Pulse Train 1';
            u_basic = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
            fcode = 0;
            F = 15;
            K = 80;
            T = 1.05;
            N = 100;
            sr = 10;

        end;


        if sigcode == 0 || sigcode == 10

            signal_name='Pulse Train 2';
            u_basic = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
            fcode = 0;
            F = 12;
            K = 80;
            T = .042;
            N = 60;
            sr = 10;

        end;

        if sigcode == 0 || sigcode == 11

            signal_name='Stepped Freq. Pulse Train';
            u_basic = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
            fcode = 1;
            f_basic = .78 * [0 0 0 0 0 1 0 0 0 0 2 0 0 0 0 3 0 0 0 0 4 0 0 0 0 5];
            F = 12;
            K = 80;
            T = .042;
            N = 60;
            sr = 10;

        end;

        if sigcode == 0 || sigcode == 12

            signal_name='Weighted Stepped Freq. Pulse Train'

            u_basic = sqrt(chebwin(36, 50))'.*[1 0 0 0 0 1 0 0 0 0 1  0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
            fcode = 1;
            f_basic = .7 * [0 0 0 0 0 1 0 0 0 0 2 0 0 0 0 3 0 0 0 0 4 0 0 0 0 5 0 0 0 0 6 0 0 0 0 7];
            F = 16;
            K = 70;
            T = .03;
            N = 60;
            sr = 5;

        end;

        ambiguity_1_scratch(u_basic,
                            fcode,
                            f_basic,
                            F,
                            K,
                            T,
                            N,
                            sr,
                            plotdir,
                            signal_name);

    end;
endfunction


