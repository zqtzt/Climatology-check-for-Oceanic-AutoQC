%%% interpolation to the standard levels, for each month
%%%%%%%%%%% 这里增加一个不插值的数组文件
clear
cd ./WOD_afterQC
file_WOD=ls;
cd ../


%Cheng的 119层
Std_depth=[0,1,2,3,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,425,450,475,500,525,550,575,600,625,650,675,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000];
Std_depth=Std_depth';
% Std_depth_viktor=[0;5;10;15;20;25;30;35;42;50;62;78;100;122;150;182;218;258;300;350;402;458;518;582;650;700;798;878;962;1050;1142;1238;1338;1442;1550;1662;1778;1898;2022;2150;2282;2418;2558;2702;2850;3002;3158;3318;3482;3650;3822;3998;4178;4362;4550;4742;4938;5138;5342;5550;5762;5978;6198];


for month=1:6
    total_number=0;
    month
    for i=3:length(file_WOD(1:end,:))
        str_cmp=['m_',num2str(month),'.mat'];
        if(contains(file_WOD(i,:),str_cmp))   %所有1月的，合并在一个数组里
            eval(['load ','./WOD_afterQC/',file_WOD(i,:),' xbtn']);
            %         eval(['load ','./WOD_afterQC/',file_WOD(i,:),' xbtn XBTdataD_new XBTdataT_new year month XBTinfo']);
            total_number=total_number+xbtn;
        end
    end
    total_number
    DataT_interp=NaN(length(Std_depth),total_number);flag
%     DataD_interp=NaN(900,total_number);
    Datainfo=NaN(18,total_number);
    
    
    %读取特定月份的数组，然后插值，最后合并
    order=1;
    for i=3:length(file_WOD(1:end,:))
        str_cmp=['m_',num2str(month),'.mat'];
        if(contains(file_WOD(i,:),str_cmp))
%             i
            eval(['load ','./WOD_afterQC/',file_WOD(i,:),' xbtn XBTdataD_new XBTdataT_new year month XBTinfo']);
            %%% 先做插值，插值到标准深度
%             T_interp(1:length(Std_depth),1:xbtn)=NaN;
            T_interp=NaN(length(Std_depth),xbtn);
            for j=1:xbtn  %循环剖面做插值
              
                depth=XBTdataD_new(:,j);
                tem=XBTdataT_new(:,j);
                depth(tem>=35 | tem <= -2.5)=NaN;
                tem(tem>=35 | tem <= -2.5)=NaN;
                tem=tem(~isnan(depth));
                depth=depth(~isnan(depth));
                ndat=length(depth);
                if(isempty(depth) | isempty(tem))
                    depth_ignore=[];
                else
                    depth_ignore=check_depth_interval(depth,Std_depth);
                end
                for std_levels=1:length(Std_depth)
                    if(~ismember(std_levels,depth_ignore))
                    T_interp(std_levels,j)=vertical_interpolation_ligc(Std_depth(std_levels),depth,tem,ndat);
                    end
                end
            end
            %%%%  然后合并数组
%             DataD_interp(:,order:order+xbtn-1)=XBTdataD_new;
            DataT_interp(:,order:order+xbtn-1)=T_interp;
            Datainfo(:,order:order+xbtn-1)=XBTinfo;
            order=order+xbtn;
            clear T_interp
        end
    end
    eval(['save ./WOD_interp/interp_m_',num2str(month),'.mat Std_depth DataT_interp Datainfo total_number -V7.3']);
    clear DataT_interp Datainfo total_number order
end

% T_interp(1:length(Std_depth),1:xbtn)=NaN;
% for i=1:xbtn
%     depth=XBTdataD_new(:,i);
%     tem=XBTdataT_new(:,i);
%     tem=tem(~isnan(depth));
%     depth=depth(~isnan(depth));
%     ndat=length(depth);
%     for std_levels=1:length(Std_depth)
%         T_interp(std_levels,i)=vertical_interpolation(Std_depth(std_levels),depth,tem,ndat);
%     end
% end

% plot(DataT_interp(1:364,1:500000),Std_depth);hold on

%%
XBTdataT(XBTdataT>35 | XBTdataT<-2)=NaN;
figure()
hold on
plot(XBTdataD(:,1:xbtn),XBTdataT(:,1:xbtn));
figure()
hold on
plot(Std_depth,DataT_interp(:,1:xbtn));