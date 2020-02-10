rngeMHz=333;
NGRD=128;
omarr=2*pi*linspace(-rngeMHz,rngeMHz,NGRD);
for B0=[34050,3300,12540,6100]
    for Rpll = 5 %6?
        strng=strcat(num2str(B0),'_',num2str(Rpll));
        for angle2="90.0" %["0.000001","45.0","90.0"]
            for cA = 0.5 %[0.1,0.5,0.9]
                for ksym = 1 %[0,0.1,1]
                    zzmax = 0;
                    for ntmix = [0.000001,100,200,500,1000]
                        a0r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(0),'_cA',num2str(cA),'_r.txt'));
                        a0i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(0),'_cA',num2str(cA),'_i.txt'));
                        a1r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                        a1i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                        s = sum(sum(a1r.*a0r+a1i.*a0i))/sum(sum(a0r.*a0r+a0i.*a0i));
                        zzmax=max(zzmax, max(max(abs(a1r-s*a0r+1i*(a1i-s*a0i)))));
                        
                    end
                    disp(zzmax)
                    for ntmix = [0.000001,100,200,500,1000]
                        a0r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(0),'_cA',num2str(cA),'_r.txt'));
                        a0i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(0),'_cA',num2str(cA),'_i.txt'));
                        a1r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                        a1i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                        s = sum(sum(a1r.*a0r+a1i.*a0i))/sum(sum(a0r.*a0r+a0i.*a0i));
                        disp(s);
                        mesh(-omarr/(2*pi),-omarr/(2*pi),abs(a1r-s*a0r+1i*(a1i-s*a0i)));%have to invert omarr
                        view([1,-1,0.5])
                        xlabel('f_1 (MHz)');
                        ylabel('f_2 (MHz)');
                        zlabel('Intensity (a.u.)')
                        zlim([0,zzmax]);
                        set(gca,'FontSize',16);
                        angle2str = num2str(round(str2double(angle2),1));
                        title(strcat('$k_{sym}=',num2str(ksym),'\times 10^6\ s^{-1},\ T_{mix}=',num2str(round(ntmix,0)),'\ ns,\ \beta_d=',angle2str,'^{\circ}$'),'Interpreter','latex','FontSize',20);
                        saveas(gcf,strcat('Chop4Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2str,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.fig'));
                        saveas(gcf,strcat('Chop4Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2str,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.pdf'));
                        close(gcf);
                    end
                end
            end
        end
    end
end