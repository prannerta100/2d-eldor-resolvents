for Rpll = 5
    for B0 = 34050
        strng = strcat(num2str(B0),'_',num2str(Rpll));
        a0r = load(strcat('Paper1noexch_spech_mat_',strng,'_Tmix1e-06ns_r.txt'));
        a0i = load(strcat('Paper1noexch_spech_mat_',strng,'_Tmix1e-06ns_i.txt'));
        a1r = load(strcat('Paper1noexch_spech_mat_',strng,'_Tmix100ns_r.txt'));
        a1i = load(strcat('Paper1noexch_spech_mat_',strng,'_Tmix100ns_i.txt'));
        s = sum(sum(a1r.*a0r+a1i.*a0i))/sum(sum(a0r.*a0r+a0i.*a0i));
        disp(s)
        disp(max(max(abs(a0r+1i*a0i))))
        mesh(abs(a1r-s*a0r+1i*(a1i-s*a0i)))
    end
end

%%
zzmax=0;
for B0=[34050,3300,12540,6100]
    for Rpll = 5 %6?
        strng=strcat(num2str(B0),'_',num2str(Rpll));
        for angle2=["0.000001","45.0","90.0"]
            for cA = [0.1,0.5,0.9]
                for ksym = [0,0.1,1]
                    a0r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(0.000001),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                    a0i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(0.000001),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                    for ntmix = [0.000001,100,200,500,1000]
                        a1r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                        a1i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                        s = sum(sum(a1r.*a0r+a1i.*a0i))/sum(sum(a0r.*a0r+a0i.*a0i));
                        zzmax=max(zzmax,max(max(abs(a1r-s*a0r+1i*(a1i-s*a0i)))));%have to invert omarr
                        
                    end
                end
            end
        end
    end
end
%%
rngeMHz=333;
NGRD=128;
omarr=2*pi*linspace(-rngeMHz,rngeMHz,NGRD);
for B0=[34050,3300,12540,6100]
    for Rpll = 5 %6?
        strng=strcat(num2str(B0),'_',num2str(Rpll));
        for angle2=["0.000001","45.0","90.0"]
            for cA = [0.1,0.5,0.9]
                for ksym = [0,0.1,1]
                    a0r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(0.000001),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                    a0i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(0.000001),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                    for ntmix = [0.000001,100,200,500,1000]
                        a1r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                        a1i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                        s = sum(sum(a1r.*a0r+a1i.*a0i))/sum(sum(a0r.*a0r+a0i.*a0i));
                        zzmax=max(max(abs(a1r+1i*a1i)));
                        disp(s);
                        mesh(-omarr/(2*pi),-omarr/(2*pi),abs(a1r-s*a0r+1i*(a1i-s*a0i)));%have to invert omarr
                        view([1,-1,0.5])
                        xlabel('f_1 (MHz)');
                        ylabel('f_2 (MHz)');
                        zlim([0,zzmax]);
                        set(gca,'FontSize',16);
                        title(strcat('$k_{sym}=',num2str(ksym*1e6),'\ s^{-1}, T_{mix}=',num2str(ntmix),' ns, \beta_d=',angle2,'$'),'Interpreter','latex','FontSize',20);
                        saveas(gcf,strcat('ChopPaper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.fig'));
                        saveas(gcf,strcat('ChopPaper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.pdf'));
                        close(gcf);
                    end
                end
            end
        end
    end
end
