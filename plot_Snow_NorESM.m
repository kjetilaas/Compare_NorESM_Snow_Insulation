%Simple plot of snow depth and soil temperature

close all
clearvars

figure; hold on;
for norESM=1:2

    if norESM==1
        filename='SnowTsoi_NorESM1_1981-2004_NH.nc';
    else
        filename='SnowTsoi_NorESM2_1980-2014_NH.nc';
    end

    snowdp=load_nc_var(filename,'SNOWDP');
    tbot=load_nc_var(filename,'TBOT');
    tsoi_full=load_nc_var(filename,'TSOI');

    if norESM==1
        layW=0.1308;  %weight for interpolation to 20cm
    else
        layW=0.6; %weight for interpolation to 20cm
    end
    tsoi=tsoi_full(:,:,1,:).*layW+tsoi_full(:,:,1,:).*(1-layW);

    %h2osno=load_nc_var(filename,'H2OSNO');
    ArrSDef_all=[]; ArrAnorm_all=[];
    for yrnr=1:size(tbot,3)/12-1
        winmonlist=(yrnr-1)*12+10:yrnr*12+3; %winter months, oct-mar
        yrmonlist=(yrnr-1)*12+10:yrnr*12+9;  %all months, oct-sep. 
        tbmean=nanmean(tbot(:,:,winmonlist),3);
        tsmean=nanmean(tsoi(:,:,winmonlist),3);
        
        %Use whole year or only winter months for amplitude calc. (probably
        %only winter used by Slater et al., 2017).
        %Aair=max(tbot(:,:,yrmonlist),[],3)-min(tbot(:,:,yrmonlist),[],3);
        %Asoi=max(tsoi(:,:,yrmonlist),[],3)-min(tsoi(:,:,yrmonlist),[],3);
        Aair=max(tbot(:,:,winmonlist),[],3)-min(tbot(:,:,winmonlist),[],3);
        Asoi=max(tsoi(:,:,winmonlist),[],3)-min(tsoi(:,:,winmonlist),[],3);    
        Anorm=(Aair-Asoi)./Aair;

        Sdepef=zeros(size(tbmean));
        for m=1:6        
            Sdepef=Sdepef+(squeeze(snowdp(:,:,winmonlist(m))).*(7-m));
        end 
        Sdepef=Sdepef./sum(1:6);

        UseData=ones(size(tbmean));

        UseData(tbmean>272.15)=0;
        UseData(tsmean>275.65)=0;
        UseData(Aair<10)=0;
        UseData(Sdepef>1.5)=0;
        UseData(Sdepef<0.01)=0;
        %Additional criteria of snow being perennial. 
        minSD=min(snowdp(:,:,yrmonlist),[],3);
        UseData(minSD>0.01)=0;
        

        ArrSDef=reshape(Sdepef',1,[]);
        ArrAnorm=reshape(Anorm',1,[]);
        ArrUseD=reshape(UseData',1,[]);

        ArrSDef=ArrSDef(ArrUseD==1);
        ArrAnorm=ArrAnorm(ArrUseD==1);
        
        ArrSDef_all=[ArrSDef_all ArrSDef];
        ArrAnorm_all=[ArrAnorm_all ArrAnorm];

        Sdepef(UseData==0)=nan;
        Anorm(UseData==0)=nan; 
        %plot(reshape(Sdepef(UseData==1)',1,[]),reshape(Anorm(UseData==1)',1,[]),'.')
        %plot(reshape(Sdepef',1,[]),reshape(Anorm',1,[]),'.')
        %figure; plot(ArrSDef,ArrAnorm,'.')
    end

    fun = @(x,ArrSDef_all)x(1)+x(2).*(1-exp(-ArrSDef_all./x(3)));

    %---plot all points, with curve fitted to all points
    % figure; plot(ArrSDef_all,ArrAnorm_all,'k.')
    % xlabel('Effective snow depth (m)'); ylabel('Anorm');
    % x0 = [0,.1,.1];
    % x = lsqcurvefit(fun,x0,ArrSDef_all,ArrAnorm_all)
    % sds = linspace(0,0.5);
    % hold on; plot(sds,fun(x,sds),'b-')
    % axis([0,1.5,0,1])

    %plot one point pr bin, with curve fitted to selected points
    binedSD=zeros(1,10); binedAn=zeros(1,10); 
    for i=1:9;
        midpoint=(i-1)*0.05+0.025;
        binedSD(i)=mean(ArrSDef_all(abs(ArrSDef_all-midpoint)<0.025));
        binedAn(i)=mean(ArrAnorm_all(abs(ArrSDef_all-midpoint)<0.025));
    %    figure; plot(ArrSDef(binedSD(i)-0.025<ArrSDef<binedSD(i)+0.025), ...
    %        ArrAnorm(binedSD(i)-0.025<ArrSDef<binedSD(i)+0.025),'.')
    end
    i=10;
    midpoint=(i-1)*0.05+0.025;
    binedSD(i)=mean(ArrSDef_all(midpoint-0.025<ArrSDef_all));
    binedAn(i)=mean(ArrAnorm_all(midpoint-0.025<ArrSDef_all));
    %figure; plot(binedSD,binedAn,'r*'); axis([0 0.5 0 1]);

    x0 = [0,.1,.1];
    x = lsqcurvefit(fun,x0,binedSD,binedAn)
    sds = linspace(0,0.5);
    %figure;plot(binedSD,binedAn,'ko',sds,fun(x,sds),'b-'); hold on;
%    plot(sds,fun(x,sds));
    axis([0,0.5,0,1])


    %---Randomly select 35 cases in each bin. and fit curve to those. 
    ntot=100;
    x_all=zeros(3,ntot);
    for nn=1:ntot;
            subsetSD=zeros(1,350); subsetAn=zeros(1,350); 
        for i=1:10
            midpoint=(i-1)*0.05+0.025; 
            if i<10
                allinbinSD=ArrSDef_all(abs(ArrSDef_all-midpoint)<0.025);
                allinbinAn=ArrAnorm_all(abs(ArrSDef_all-midpoint)<0.025);
            else
                allinbinSD=ArrSDef_all(midpoint-0.025<ArrSDef_all);
                allinbinAn=ArrAnorm_all(midpoint-0.025<ArrSDef_all);
            end
            for ii=(i-1)*35+1:i*35
                iirand=randi(length(allinbinSD));
                subsetSD(ii)=allinbinSD(iirand);
                subsetAn(ii)=allinbinAn(iirand);
            end
        end
        x0 = [0,.1,.1];
        x = lsqcurvefit(fun,x0,subsetSD,subsetAn)
        x_all(:,nn)=x;
        sds = linspace(0,0.5);
        %plot(subsetSD,subsetAn,'.',sds,fun(x,sds),'--')
    end
    xmedian=median(x_all,2);
    plot(sds,fun(xmedian,sds),'-','LineWidth',2)
    axis([0,0.5,0,1])
end %NorESM

%legend('NorESM1 all','NorESM1 subset','NorESM2 all','NorESM2 subset');
legend('NorESM1','NorESM2');
xlabel('Effective snow depth (m)'); ylabel('Anorm');
saveas(gcf, 'Snow_Insulation.png');