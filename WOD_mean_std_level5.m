%%%%���������ϵ��ֳ��۲����ݣ����ָ�㣬Ȧ��ÿһ������еĹ۲����ݣ�Ȼ����ÿһ��ı�׼�ƽ��ֵʲô�� ������1x1���reference����
%%% mean_std_level
clear
clc

% ��1x1�ȵ����񣬹�������̬�����ĸ��Ϊ����λ��������Ϊ0.5�ȵ�С��λ��
lat=-89:90;
lon=-179:180;
lat_bnd=[lat-0.5; lat+0.5];
lon_bnd=[lon-0.5; lon+0.5];

load('F:\Matlabѧϰ\mask_basin_1_new.mat')
%%%%��ȡ�۲�����
for month=7:12
    month
    clear T_stastical*
    clear Std_depth
    %1��
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
            %%%����Ҫ�޸ģ������һ����㣬��׼���쳣�ô󣬶����������٣������ð뾶����������ĸ������ݣ���ô����϶��ǰ����˲�ͬ��ˮ�ţ��Ǵ���ģ�����Ҫ�����ģ������ķ����ǲ������������ˣ�������ֱ�Ӻ��Ե���������û�й۲�
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
            number_obs=histc(double(index),1);  %��������
            number_obs=number_obs-histc(single(all(isnan(temp))==1),1);  %��ȫΪNaN������������޳���
            if(number_obs==0 && mask(i,j)==0)  %�����һ���޹۲⣬����½���ϣ��Ͳ����ˣ�ֱ������
                %             T_stastical_mean(i,j,1:63)=NaN;
                %             T_stastical_std(i,j,1:63)=NaN;
                %             T_stastical_number(i,j,1:63)=NaN;
                continue
            end
            
            %%%����Ķ����ں��ϵĵ㣬�����ǽ����ĵ�
            % �����޹۲�ĵ㣬�ͻ�ԲȦ��������
            
            %%%%ÿһ��һ��ƽ��ֵ��stdֵ
            if(number_obs>=5)  %�����ﵽ5������ֱ�Ӽ���   %����1000m���¿��ܻ�������
                c=~isnan(temp);
                number=histc(double(c)',1)';
                [temp,number]=remove_max_min(number,temp,8);  %����8���۲�ֵ��level��ȥ��һ�����ֵһ����Сֵ
                T_stastical_number(i,j,:)=number;
                T_stastical_mean(i,j,:)=nanmean(temp');
                T_stastical_std(i,j,:)=nanstd(temp');
                %                 figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=0 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
            else
                %%%%%%%%����Ҫ���������������������ٵģ������ﲻ��5���ģ�����Χ�������ݣ�1���ڵĻ���2���ڵģ�������
                distance=sqrt((Datainfo(6,:)-lat(j)).^2+(Datainfo(7,:)-lon(i)).^2);
                index2=distance<=1;
                temp=DataT_interp(:,index2);
                number_obs=histc(double(index2),1);  %��������
                flag=1;
                if(number_obs<5)
                    %�뾶����1.5��
                    distance=sqrt((Datainfo(6,:)-lat(j)).^2+(Datainfo(7,:)-lon(i)).^2);
                    index2=distance<=1.5;
                    temp=DataT_interp(:,index2);
                    number_obs=histc(double(index2),1);  %��������
                    
                    if(number_obs <5)  %�뾶����1.5��
                        index2=distance<=2;
                        temp=DataT_interp(:,index2);
                        number_obs=histc(double(index2),1);  %��������
                        
                        if(number_obs==0) %Ҫ�ǻ�������5���Ļ�����׼��ֱ��Ϊ0���ͷ���������   %%�������˼����д����
                            T_stastical_number(i,j,:)=0;
                            continue
                        elseif(number_obs==1)
                            T_stastical_number(i,j,:)=double(~isnan(temp));
                            T_stastical_mean(i,j,:)=temp;
                            temp_std=temp;
                            temp_std(:)=0;
                            T_stastical_std(i,j,:)=temp_std; %ֻ��һ��������׼��Ϊ0
                            continue
                        end
                        c=~isnan(temp);
                        number=histc(double(c)',1)';
                        [temp,number]=remove_max_min(number,temp,8);  %����8���۲�ֵ��level��ȥ��һ�����ֵһ����Сֵ
                        T_stastical_number(i,j,:)=number;
                        T_stastical_mean(i,j,:)=nanmean(temp');
                        T_stastical_std(i,j,:)=nanstd(temp');
                        %                     figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=3 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
                    end
                    %��Щ���ٶ��������۲�������
                    c=~isnan(temp);
                    number=histc(double(c)',1)';
                    [temp,number]=remove_max_min(number,temp,8);  %����8���۲�ֵ��level��ȥ��һ�����ֵһ����Сֵ
                    T_stastical_number(i,j,:)=number;
                    T_stastical_mean(i,j,:)=nanmean(temp');
                    T_stastical_std(i,j,:)=nanstd(temp');
                    %                     figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=2 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
                    continue
                end
                c=~isnan(temp);
                number=histc(double(c)',1)';
                [temp,number]=remove_max_min(number,temp,8);  %����8���۲�ֵ��level��ȥ��һ�����ֵһ����Сֵ
                T_stastical_number(i,j,:)=number;
                T_stastical_mean(i,j,:)=nanmean(temp');
                T_stastical_std(i,j,:)=nanstd(temp');
                %                 figure();plot(temp,Std_depth,'-');set(gca,'YDir','reverse');title(['flag=1 lon=',num2str(lon(i)),' lat=',num2str(lat(j))]);pause;close(gcf);
            end
            
        end
    end
    
    eval(['save ./WOD_stastical_field/WOD_stastical_mean_std_number_month',num2str(month),'_new.mat lat lon lat_bnd lon_bnd T_stastical* Std_depth']);
    
    %%%% ��ֵ�˲� 9��ռ��˲�
    for k=1:length(Std_depth)
        T_stastical_mean(:,:,k)=movmedian2d(movmedian2d(T_stastical_mean(:,:,k),9),9);
        T_stastical_std(:,:,k)=movmedian2d(movmedian2d(T_stastical_std(:,:,k),9),9);
    end
    eval(['save ./WOD_stastical_field/after_medfilt/WOD_stastical_mean_std_number_month',num2str(month),'.mat lat lon lat_bnd lon_bnd T_stastical* Std_depth']);
end
a='ok'
%%%%��ȫ���λ ƽ��ֵ����׼��ֲ�ͼ

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
%%%%%%%%%%%%  ����Χ9���ֵ����亣��ΪNaN�ĵ㣬 ����Ӧ�ĳɲ�ֵ�����������·ݵ����ݣ�
load mask_basin_1.mat
for k=1:38 %ֻ���ϲ�2000m�ķ�½�صĵ����9��ƽ��
    k
    T_mean_level=T_stastical_mean(:,:,k);
    T_std_level=T_stastical_std(:,:,k);
    [row,col]=find(mask~=0);
    for i=1:length(row)
        if(isnan(T_mean_level(row(i),col(i))))  %������ܳ���
            %����� row(i),col(i)
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
m_grid('linestyle','none','box','tickdir','out')  %����γ������
m_coast('patch',0.95*[207 207 207]/255,'edgecolor',0.95*[207 207 207]/255);
caxis([-2 35])
colorbar
load red_blue_16.mat
% load New_England.mat
% title('300m')
colormap(mycolor)