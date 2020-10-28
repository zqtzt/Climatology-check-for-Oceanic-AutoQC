%%% 用WOA18的方法，在5x5读的格点，用莱茵达准则，去掉3~5倍标准差
clear
%先建立5x5°的格点，中心为2.5度
lat_5=-87.5:5:87.5;
lat_5_bnd=[lat_5-2.5;lat_5+2.5];
lon_5=-177.5:5:177.5;
lon_5_bnd=[lon_5-2.5;lon_5+2.5];

load ./WOD_Rhenida/mask_5.mat  %1包含陆地，2靠近陆地，3开放大洋+陆地


for month=2:12
    month
    eval(['load ./WOD_interp/未经过Rehinda准则但是经过表层平均的/interp_m_',num2str(month),'.mat']);
    
    T_number=NaN(72,36,120);
    T_mean=NaN(72,36,120);
    T_std=NaN(72,36,120);
    
    
    
    %每一个格点都统计出来
    for i=1:length(lon_5)
        i
        for j=1:length(lat_5)
            flag=0;
            index=(Datainfo(6,:) >= lat_5_bnd(1,j) & Datainfo(6,:) <= lat_5_bnd(2,j)) & (Datainfo(7,:) >= lon_5_bnd(1,i) & Datainfo(7,:) <= lon_5_bnd(2,i));
            temp=DataT_interp(:,index);
            temp(temp>35 | temp <-2)=NaN;
            c=~isnan(temp);
            number=histc(double(c)',1)';
            T_number(i,j,:)=number;
            T_mean(i,j,:)=nanmean(temp');
            T_std(i,j,:)=nanstd(temp');
        end
    end
    
    %这个用来标记50-400m，变率大的格点
    flag_large_variablity=NaN(72,36,120);
    for k=1:length(Std_depth)
        if(k<=12 || k >=100)
            Q80_std(k)=NaN;
        else
            a=T_std(:,:,k);
            Q80_std(k)=prctile(a(:),85);
            flag_large_variablity(:,:,k)=single(a>Q80_std(k));
        end
    end
    
    eval(['save ./WOD_Rhenida/T_stat_5_month',num2str(month),'.mat T_number T_mean T_std lon_5 lat_5 Std_depth lat_5_bnd lon_5_bnd'])
    
    for m=1:length(DataT_interp)
%         m
        %     find(Datainfo(6,m) >= lat_5_bnd(1,:) & Datainfo(6,m) <= lat_5_bnd(2,:)) & (Datainfo(7,m) >= lon_5_bnd(1,:) & Datainfo(7,m) <= lon_5_bnd(2,:))
        lat_index=find(Datainfo(6,m) >= lat_5_bnd(1,:) & Datainfo(6,m) < lat_5_bnd(2,:));
        lon_index=find(Datainfo(7,m) >= lon_5_bnd(1,:) & Datainfo(7,m) < lon_5_bnd(2,:));
        try   %对于上层50m
            switch mask_5(lon_index,lat_index)
                case 3  %开放大洋 >std
                    thresold_right=reshape(T_mean(lon_index,lat_index,1:12)+3*T_std(lon_index,lat_index,1:12),12,1);
                    thresold_left=reshape(T_mean(lon_index,lat_index,1:12)-3*T_std(lon_index,lat_index,1:12),12,1);
                    DataT_interp(DataT_interp(1:12,m) > thresold_right | DataT_interp(1:12,m) < thresold_left,m)=NaN;
                case 2 %近海海岸 >4std
                    thresold_right=reshape(T_mean(lon_index,lat_index,1:12)+4*T_std(lon_index,lat_index,1:12),12,1);
                    thresold_left=reshape(T_mean(lon_index,lat_index,1:12)-4*T_std(lon_index,lat_index,1:12),12,1);
                    DataT_interp(DataT_interp(1:12,m) > thresold_right | DataT_interp(1:12,m) < thresold_left,m)=NaN;
                case 1  %包含陆地格点， >5std
                    thresold_right=reshape(T_mean(lon_index,lat_index,1:12)+5*T_std(lon_index,lat_index,1:12),12,1);
                    thresold_left=reshape(T_mean(lon_index,lat_index,1:12)-5*T_std(lon_index,lat_index,1:12),12,1);
                    DataT_interp(DataT_interp(1:12,m) > thresold_right | DataT_interp(1:12,m) < thresold_left,m)=NaN;
            end
        catch
            m,1
        end
        %对于50m-4000m
        try
            for k=13:99
                switch flag_large_variablity(lon_index,lat_index,k)
                    case 1 %变率较大（大于85%分位数）的格点
                        thresold_right=T_mean(lon_index,lat_index,k)+5*T_std(lon_index,lat_index,k);
                        thresold_left=T_mean(lon_index,lat_index,k)-5*T_std(lon_index,lat_index,k);
                        if (DataT_interp(k,m) > thresold_right || DataT_interp(k,m) < thresold_left)
                            DataT_interp(k,m)=NaN;
                        end
                    case 0 %变率较小的格点
                        thresold_right=T_mean(lon_index,lat_index,k)+3*T_std(lon_index,lat_index,k);
                        thresold_left=T_mean(lon_index,lat_index,k)-3*T_std(lon_index,lat_index,k);
                        if (DataT_interp(k,m) > thresold_right || DataT_interp(k,m) < thresold_left)
                            DataT_interp(k,m)=NaN;
                        end
                end
            end
        catch
            hh='I am here'
        end
    end
    
    eval(['save ./WOD_interp/interp_m_',num2str(month),'.mat Datainfo DataT_interp Std_depth total_number'])
end
% hold on
% [lon1,lat1]=meshgrid(lon_5-2.5,lat_5-2.5);
% m_proj('miller','lon',[-180 180],'lat',[-90 90]);
% order=0;
% for j=1:length(lat_5)
%     for i=1:length(lon_5)
%         order=order+1;
%         %         m_plot(lon_5(i),lat_5(j),'Marker','o','Color','black','LineStyle','none','MarkerSize',5,'MarkerFaceColor','black')
%         m_text(lon_5(i),lat_5(j),num2str(order),'HorizontalAlignment','center','FontSize',8)
%         m_rectangle(lon_5(i)-2.5,lat_5(j)-2.5,5,5,'Curvature',[0,0]);
%     end
% end
% m_pcolor(lon1,lat1,mask_5','linestyle','none');
% m_grid('box','fancy','linestyle','none','gridcolor','w','backcolor',[0.2 0.65 1]);
% m_coast('patch',[0.85 0.33 0.1],'edgecolor','none');

