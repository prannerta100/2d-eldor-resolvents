rngeMHz=333;
NGRD=128;
omarr=2*pi*linspace(-rngeMHz,rngeMHz,NGRD);
for B0=[34050,3300,12540,6100]
    for Rpll = 5 %6?
        strng=strcat(num2str(B0),'_',num2str(Rpll));
        for angle2=["0.000001","45.0","90.0"]
            for cA = [0.1,0.5,0.9]
                for ksym = [0,0.1,1]
                    for ntmix = [500]%[0.000001,100,200,500,1000]
                        a1r=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_r.txt'));
                        a1i=load(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'_i.txt'));
                        plot(-omarr(1:end-1)/(2*pi),diff(sum(a1r,2)));
                        hold on;
                        plot(-omarr(1:end-1)/(2*pi),diff(sum(a1r,1)));
                        hold off;
                        xlabel('f_1 (MHz)');
                        ylabel('f_2 (MHz)');
                        set(gca,'FontSize',16);
                        angle2str = num2str(round(str2double(angle2),1));
                        title(strcat('$k_{sym}=',num2str(ksym),'\times 10^6\ s^{-1}, T_{mix}=',num2str(ntmix),' ns, \beta_d=',angle2str,'$'),'Interpreter','latex','FontSize',20);
                        saveas(gcf,strcat('Chop3Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2str,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.fig'));
                        saveas(gcf,strcat('Chop3Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',angle2str,'deg_ksym',num2str(ksym*1e6),'_cA',num2str(cA),'.pdf'));
                        close(gcf);
                    end
                end 
            end
        end
    end
end