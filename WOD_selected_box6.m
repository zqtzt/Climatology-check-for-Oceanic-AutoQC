clear

for month=7:12
    month
    clear t_max t_min flag_row* flag_column* T_stastical*
    eval(['load F:\QC_science\climatology_study\WOD_observations_2020up\WOD_stastical_field\after_medfilt\WOD_stastical_mean_std_number_month',num2str(month),'.mat']);
    
    %%%% 读取每个格点的平均态，按照标准差选格点
    t_max=T_stastical_mean+1.*T_stastical_std;
    t_min=T_stastical_mean-1.*T_stastical_std;
    
    flag_row=cell(360,180,79);
    flag_column=cell(360,180,79);
    
    flag_row_inversion=cell(360,180,79);  %记录反选的格点
    flag_column_inversion=cell(360,180,79);
    
    flag_raidus=NaN(360,180,79);  %用于标注中心格点的选格点的半径 （我们只考虑上层2000m的变率较大的地方用3°）
    clear Q90
    %%%%% 上层1000m，找标准差的前90%的分位数
    for k=1:length(Std_depth)
        if(k<=60)  %这里改！！！！！！！！！！
            t_std_level=T_stastical_std(:,:,k);
            Q90(k)=prctile(t_std_level(:),90);
        else
            Q90(k)=999;
        end
    end
    
    for i=1:360 %longitude
        i
        for j=1:180 %latitude
            for k=1:79 %先做上层1950m的
                %             clear Q90
                t_std_level=T_stastical_std(:,:,k);
                t_max_level=t_max(:,:,k);
                t_min_level=t_min(:,:,k);
                flag=NaN(360,180);
                if(isnan(T_stastical_mean(i,j,k)))  %陆地的去掉  %这里要改，因为这里未必是陆地，可能是中间没有具体数值的格点！！！！
                    %flag_row{i,j,k}=NaN;
                    %flag_column{i,j,k}=NaN;
                    continue
                end
                
                %             t_std_level(i,j)
                if(t_std_level(i,j)>Q90(k))  %变率较大的地方用3°半径来选水团相似的格点
                    radius=3;
                    [flag]=select_box_condition(flag,i,j,k,radius,lat,lon,t_min_level,t_max_level,Q90,t_std_level);
                    %               %%%%% 重新建立一个数组，标注哪些格点是半径为3的，标注的在后面选数据的时候，要选临近三个月的数据
                    flag_raidus(i,j,k)=3;
                    kk=flag;
                else
                    radius=5;
                    flag_raidus(i,j,k)=5;
                    [flag]=select_box_condition(flag,i,j,k,radius,lat,lon,t_min_level,t_max_level,Q90,t_std_level);
                    
                    
                end
                [row,col]=find(flag==1);  %被选上的相邻格点
                flag_row{i,j,k}=row;  %元胞数组  行=经度i，列=维度j
                flag_column{i,j,k}=col;
                
                [row_inversion,col_inversion]=find(flag==0);
                flag_row_inversion{i,j,k}=row_inversion;  %元胞数组  行=经度i，列=维度j
                flag_column_inversion{i,j,k}=col_inversion;
                
            end
        end
    end
    
    filenames=['./WOD_选格点/upper1950/flag_infos_month',num2str(month),'_new.mat'];
    eval(['save ',filenames,' flag_column flag_row lat lon lat_bnd lon_bnd Std_depth flag_row_inversion flag_column_inversion flag_raidus']);
    % save ./WOD_选格点/flag_use_threemonths_info_month1.mat flag_raidus
    clear flag_infos
    
    
end

%%  按照全部距离，圈符合5度的格点
d='I am here'
pause
clear
load('./WOD_stastical_field/WOD_stastical_mean_std_number_month1.mat')
%%%%% 上层1000m，找标准差的前90%的分位数
for k=1:length(Std_depth)
    if(k<=30)
        t_std_level=T_stastical_std(:,:,k);
        Q90(k)=prctile(t_std_level(:),90);
    else
        Q90(k)=999;
    end
end

t_max=T_stastical_mean+1.*T_stastical_std;

flag_row_noncondition=cell(360,180,63);
flag_column_noncondition=cell(360,180,63);
load mask_basin_1.mat
for i=1:360 %longitude
    i
    for j=1:180 %latitude
        for k=2:2 %38 length(Std_depth)
            flag=NaN(360,180);
            t_std_level=T_stastical_std(:,:,k);
            if(isnan(t_max(i,j,k)))  %陆地的去掉  %这里要改，因为这里未必是陆地，可能是中间没有具体数值的格点！！！！
                %这里还要改，在下层2000m以下就不起作用了
                %flag_row{i,j,k}=NaN;
                %flag_column{i,j,k}=NaN;
                continue
            end
            if(t_std_level(i,j)>Q90(k))
                radius=3;
                [flag]=select_box_noncondition(flag,i,j,radius,lat,lon);
                kk=flag;
                %                 pause
            else
                radius=5;
                [flag]=select_box_noncondition(flag,i,j,radius,lat,lon);
            end
            
            [row,col]=find(flag==1);
            flag_row_noncondition{i,j,k}=row;  %元胞数组  行=经度i，列=维度j
            flag_column_noncondition{i,j,k}=col;
            
        end
    end
end
save ./WOD_选格点/flag_infos_month1_noncondition_boxes.mat flag_column_noncondition flag_row_noncondition lat lon lat_bnd lon_bnd


%%
clear
load ./WOD_选格点/flag_infos_month1_new.mat
load('F:\QC_science\climatology_study\WOD_observations_2020up\WOD_stastical_field\after_medfilt\WOD_stastical_mean_std_number_month1_new.mat')
% load('F:\QC_science\climatology_study\WOD_observations_2020up\WOD_stastical_field\WOD_stastical_mean_std_number_month1_new.mat')
% load ./WOD_选格点/flag_infos_month1_noncondition_boxes.mat
%%  画圆圈，画格点，看看每一个格点的情况！！！
load mycolor_England_1.mat
% mycolor=mycolor(13:end,:);
selected_grid=NaN(360,180);

%i=147;
% j=91;
% k=2;
for i=303:303
    i
    for j=120:120
        for k=3:3
            %         k=3;
            radius=flag_raidus(i,j,k);
            if(isempty(flag_row{i,j,k}))
                continue
            end
            selected_grid=NaN(360,180);
            row=flag_row{i,j,k}; %被标记的行号 列好
            column=flag_column{i,j,k};
            row_inversion=flag_row_inversion{i,j,k};
            column_inversion=flag_column_inversion{i,j,k};
            
            for m=1:length(column)
                selected_grid(row(m),column(m))=1;
            end
            
            figure()
            set(gcf,'Position',[1000,740,980,598]);
            subplot(1,2,1)
            %         set(gcf,'visible','off');
            hold on
            [lon1,lat1]=meshgrid(lon-0.5,lat-0.5);
            % m_proj('robinson','clongitude',90);
            lon2=double(round(lon(i)));
            lat2=double(round(lat(j)));
            if (lat2-10 <= -90)
                lat2=-79.5;
            elseif(lat2+10 >= 90)
                lat2=79.5;
            end
            m_proj('miller','lon',[lon2-10 lon2+10],'lat',[lat2-10 lat2+10]);
            m_pcolor(lon1,lat1,T_stastical_mean(:,:,k)','linestyle','none');
            %         gca1=m_pcolor(lon1,lat1,selected_grid','linestyle','none');
            m_plot(lon(row),lat(column),'Marker','o','Color','black','LineStyle','none','MarkerSize',5,'MarkerFaceColor','black')
            m_plot(lon(i),lat(j),'Marker','p','Color','red','LineStyle','none','MarkerSize',5,'MarkerFaceColor','red')
            m_rectangle(lon(i)-radius,lat(j)-radius,radius*2,radius*2,'Curvature',[1,1]);
            % m_rectangle(lon(i)-3,lat(j)-3,6,6,'Curvature',[1,1]);
            m_grid('box','fancy','linestyle','none','gridcolor','w','backcolor',[0.2 0.65 1]);
            m_coast('patch',[0.85 0.33 0.1],'edgecolor','none');
            total_number=histc(selected_grid(:),1);
            txt=sprintf('total grid number= %d',total_number);
            title(txt)
            %         caxis([20 33])
            colormap(mycolor)
            colorbar
            
            subplot(1,2,2)
            try
                plot_selected_boxed_range(i,j,k,row,column,row_inversion,column_inversion,T_stastical_mean,T_stastical_std)
            end
            pause
            close(gcf)
        end
    end
end
% save example_1std.mat info lat lon lat_bnd lon_bnd




%%
for n=1:length(row_inversion)
    T_mean(n)=T_stastical_mean(row_inversion(n),column_inversion(n),k);
    T_std(n)=T_stastical_std(row_inversion(n),column_inversion(n),k);
end
figure()
histogram(T_mean,10)

radius=flag_raidus(:,:,3);
radius(radius~=3)=NaN;
plot_lat_lon(lon,lat,radius,'STD>Q90, radius=3');