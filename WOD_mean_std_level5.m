%%%%利用我手上的现场观测数据，划分格点，圈出每一个格点有的观测数据，然后算每一层的标准差、平均值什么的 （构建1x1°的reference场）
%%% mean_std_level
clear
clc

% 用1x1度的网格，构建气候态，中心格点为整数位数，两侧为0.5度的小数位数
lat=-89:90;
lon=-179:180;
lat_bnd=[lat-0.5; lat+0.5];
lon_bnd=[lon-0.5; lon+0.5];

load('F:\Matlab学习\mask_basin_1_new.mat')
%%%%读取观测数据
for month=7:12
    month
    clear T_stastical*
    clear Std_depth
    %1月
    eval(['load ./WOD_interp/interp_m_',num2str(month),'.mat']);
    
    T_stastical_mean=NaN(360,180,length(Std_depth));
    T_stastical_std=NaN(360,180,length(Std_depth));
    T_stastical_number=NaN(360,180,length(Std_depth));
    
    for i=1:length(lon) %longitude
        i
        flag=0;
        for j=1:length(lat) %latitude
            %           for m=1:length(row)
            %               i=row(m);
            %               j=column(m);
            clear index temp c
            %%%这里要修改，如果上一个格点，标准差异常得大，而且数量又少，又是用半径囊括了外面的格点的数据，那么这个肯定是包括了不同的水团，是错误的，是需要订正的，订正的方法是不考虑这个格点了，这个格点直接忽略掉，本来就没有观测
            for p=1:length(Std_depth)
                if(flag==1 & T_stastical_std(i,j-1,p)>=3 & T_stastical_number(i,j-1,p)<=10)
                    T_stastical_std(i,j-1,p)=0;
                    T_stastical_mean(i,j-1,p)=NaN;
                end
            end
            flag=0;
            index=(Datainfo(6,:) >= lat_bnd(1,j) & Datainfo(6,:) <= lat_bnd(2,j)) & (Datainfo(7,:) >= lon_bnd(1,i) & Datainfo(7,:) <= lon_bnd(2,i));
            temp=DataT_interp(:,index);
            temp(temp>35 | temp <-2)=NaN;
            number_obs=histc(double(index),1);  %剖面数量
            number_obs=number_obs-histc(single(all(isnan(temp))==1),1);  %把全为NaN的剖面的数量剔除掉
            if(number_obs==0 && mask(i,j)==0)  %如果这一点无观测，且在陆地上，就不管了，直接跳过
                %             T_stastical_mean(i,j,1:63)=NaN;
                %             T_stastical_std(i,j,1:63)=NaN;
                %             T_stastical_number(i,j,1:63)=NaN;
                continue
            end
            
            %%%下面的都是在海上的点，或者是近海的点
            % 海上无观测的点，就画圆圈囊括进来
            
            %%%%每一层一个平均值和std值
            if(number_obs>=5)  %数量达到5个，则直接计算   %这在1000m以下可能还有困难
                c=~isnan(temp);
                number=histc(double(c)',1)';
                [temp,number]=remove_max_min(number,temp,8);  %大于8个观测值的level，去掉一个最大值一个最小值
                T_stastical_number(i,j,:)=number;
                T_stastical_mean(i,j,:)=nanmean(temp');
                T_stastical_std(i,j,:)=nanstd(temp');
                %                 figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=0 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
            else
                %%%%%%%%这里要加条件限制数量，数量少的，数量达不到5个的，用周围格点的数据（1度内的或者2度内的）来补充
                distance=sqrt((Datainfo(6,:)-lat(j)).^2+(Datainfo(7,:)-lon(i)).^2);
                index2=distance<=1;
                temp=DataT_interp(:,index2);
                number_obs=histc(double(index2),1);  %剖面数量
                flag=1;
                if(number_obs<5)
                    %半径换成1.5度
                    distance=sqrt((Datainfo(6,:)-lat(j)).^2+(Datainfo(7,:)-lon(i)).^2);
                    index2=distance<=1.5;
                    temp=DataT_interp(:,index2);
                    number_obs=histc(double(index2),1);  %剖面数量
                    
                    if(number_obs <5)  %半径换成1.5度
                        index2=distance<=2;
                        temp=DataT_interp(:,index2);
                        number_obs=histc(double(index2),1);  %剖面数量
                        
                        if(number_obs==0) %要是还是少于5个的话，标准差直接为0，就放弃治疗了   %%这里谨慎思考改写！！
                            T_stastical_number(i,j,:)=0;
                            continue
                        elseif(number_obs==1)
                            T_stastical_number(i,j,:)=double(~isnan(temp));
                            T_stastical_mean(i,j,:)=temp;
                            temp_std=temp;
                            temp_std(:)=0;
                            T_stastical_std(i,j,:)=temp_std; %只有一个数，标准差为0
                            continue
                        end
                        c=~isnan(temp);
                        number=histc(double(c)',1)';
                        [temp,number]=remove_max_min(number,temp,8);  %大于8个观测值的level，去掉一个最大值一个最小值
                        T_stastical_number(i,j,:)=number;
                        T_stastical_mean(i,j,:)=nanmean(temp');
                        T_stastical_std(i,j,:)=nanstd(temp');
                        %                     figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=3 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
                    end
                    %这些至少都是两个观测剖面了
                    c=~isnan(temp);
                    number=histc(double(c)',1)';
                    [temp,number]=remove_max_min(number,temp,8);  %大于8个观测值的level，去掉一个最大值一个最小值
                    T_stastical_number(i,j,:)=number;
                    T_stastical_mean(i,j,:)=nanmean(temp');
                    T_stastical_std(i,j,:)=nanstd(temp');
                    %                     figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=2 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
                    continue
                end
                c=~isnan(temp);
                number=histc(double(c)',1)';
                [temp,number]=remove_max_min(number,temp,8);  %大于8个观测值的level，去掉一个最大值一个最小值
                T_stastical_number(i,j,:)=number;
                T_stastical_mean(i,j,:)=nanmean(temp');
                T_stastical_std(i,j,:)=nanstd(temp');
                %                 figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=1 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
            end
            
        end
    end
    
    eval(['save ./WOD_stastical_field/WOD_stastical_mean_std_number_month',num2str(month),'_new.mat lat lon lat_bnd lon_bnd T_stastical* Std_depth']);
    
    %%%% 中值滤波 9点空间滤波
    for k=1:length(Std_depth)
        T_stastical_mean(:,:,k)=movmedian2d(movmedian2d(T_stastical_mean(:,:,k),9),9);
        T_stastical_std(:,:,k)=movmedian2d(movmedian2d(T_stastical_std(:,:,k),9),9);
    end
    eval(['save ./WOD_stastical_field/after_medfilt/WOD_stastical_mean_std_number_month',num2str(month),'.mat lat lon lat_bnd lon_bnd T_stastical* Std_depth']);
end
a='ok'
%%%%画全球二位 平均值、标准差分布图

% plot_lat_lon(lon,lat,result2,'Reference field T\_stastical\_std');
% load rainbow.mat
% load red_blue_16.mat
% colormap(mycolor)
% caxis([0 6])
% % %
% % for k=99:length(Std_depth)
% plot_lat_lon(lon,lat,result,'Reference field T\_stastical\_mean');
% set(gcf,'Position',[1000,605,944,733])
% load rainbow.mat
% load red_blue_16.mat
% colormap(mycolor)
% caxis([-2 35])
% pause
% close(gcf)
% end
% %
% %
% c=T_stastical_number(:,:,7);
% c(c>=3)=NaN;
% plot_lat_lon(lon,lat,c,'Reference field T\_stastical\_number');
% load rainbow.mat
% % load red_blue_16.mat
% colormap(mycolor(:,1:3))
% caxis([0 2])


% global_weight_ave_level(T_stastical_std(:,:,3),lat,lon)
%%
%%%%%%%%%%%%  用周围9点均值，填充海上为NaN的点， 后面应改成插值，或者相邻月份的数据！
load mask_basin_1.mat
for k=1:38 %只对上层2000m的非陆地的点进行9点平均
    k
    T_mean_level=T_stastical_mean(:,:,k);
    T_std_level=T_stastical_std(:,:,k);
    [row,col]=find(mask~=0);
    for i=1:length(row)
        if(isnan(T_mean_level(row(i),col(i))))  %这里可能出错
            %坐标点 row(i),col(i)
            try
                temp=[T_mean_level(row(i)-1,col(i)+1) T_mean_level(row(i)-1,col(i)) T_mean_level(row(i)-1,col(i)-1) T_mean_level(row(i),col(i)-1) ...
                    T_mean_level(row(i),col(i)+1) T_mean_level(row(i)+1,col(i)+1) T_mean_level(row(i)+1,col(i)) T_mean_level(row(i)+1,col(i)-1)];
                temp_std=[T_std_level(row(i)-1,col(i)+1) T_std_level(row(i)-1,col(i)) T_std_level(row(i)-1,col(i)-1) T_std_level(row(i),col(i)-1) ...
                    T_std_level(row(i),col(i)+1) T_std_level(row(i)+1,col(i)+1) T_std_level(row(i)+1,col(i)) T_std_level(row(i)+1,col(i)-1)];
            catch
                continue
            end
            T_stastical_mean(row(i),col(i),k)=nanmean(temp);
            T_stastical_std(row(i),col(i),k)=nanmean(temp_std);
            
        end
    end
end

save ./WOD_stastical_field/WOD_stastical_mean_std_number_month1.mat lat lon lat_bnd lon_bnd T_stastical* Std_depth



%%
load ('./WOD_stastical_field/WOD_stastical_mean_std_number_month1.mat')
filename='../WOA18_1/month/woa18_decav_t01_01.nc';
file_keys = netcdf.open(filename,'NOWRITE');
[t_number,~]=get_nc_variable(file_keys,'t_dd');
[t_mean,fill_value]=get_nc_variable(file_keys,'t_mn');
[t_std,fill_value]=get_nc_variable(file_keys,'t_sd');
[depth_WOA,~]=get_nc_variable(file_keys,'depth');

%%
[lon1,lat1]=meshgrid(double(lon),double(lat));
% t_number=double(t_number);
% T_diff=T_stastical_number(:,:,2)-t_number(:,:,2);

figure()
set(gcf,'Position',[488,219.4,762.5999999999999,542.6]);
% m_proj('robinson','clongitude',180);
m_proj('robinson','lon',[-180 180],'lat',[-90 90]);
gca1=m_pcolor(lon1,lat1,T_stastical_mean(:,:,1)','linestyle','none');
%  gca1=m_pcolor(lon1,lat1,t_number(:,:,2)','linestyle','none');
m_grid('linestyle','none','box','tickdir','out')  %画经纬度坐标
m_coast('patch',0.95*[207 207 207]/255,'edgecolor',0.95*[207 207 207]/255);
caxis([-2 35])
colorbar
load red_blue_16.mat
% load New_England.mat
% title('300m')
colormap(mycolor)