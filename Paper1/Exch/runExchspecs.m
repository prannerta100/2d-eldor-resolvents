for B0=34050%[3300,34050,12540,6100]
    for Rpll = [5,6]
        for angle2=["45.0","90.0"]
            for cA = 0.5
                for ksym = [0,1,10]
                    for ntmix = [0.000001,100,200,500,1000]
                        spec_mat=expmeldor_exchange_prune(cA,ksym,Rpll,B0,ntmix,'0.000001',angle2);
                    end
                end
            end
        end
    end
end
%%
for B0=[34050,3300,12540,6100]
    for Rpll = 5
        for angle2=["0.000001","45.0","90.0"]
            for cA = [0.1,0.5,0.9]
                for ksym = [0,0.1,1]
                    ntmixlist = [0.000001,100,200,500,1000];
                    spec_mat=expmeldor_exchange_prune_fast(cA,ksym,Rpll,B0,ntmixlist,'0.000001',angle2);
                end
            end
        end
    end
end