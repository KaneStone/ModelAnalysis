%%testing why m_map hatching doesnt work

testdata = rand(145,97);
longitude = 0:2.5:360;
latitude = -90:1.875:90;

m_proj('Stereographic','lon',0,'lat',90,'rad',abs(90))                

[~,h] = m_contourf(longitude,latitude,testdata','LineStyle','none');        

m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);


bndry_lon=[-170 -120 -120 -170 -170];
bndry_lat=[49 49 70 70 49];

%m_line(bndry_lon,bndry_lat,'linewi',2,'color','k');     % Area outline ...
m_hatch(bndry_lon,bndry_lat,'single',1,5,'color','k'); % ...with hatching added.