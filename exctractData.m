
    clc
    figName = '160328-IMG-a-fit.fig';
    [pathstr,name,ext] = fileparts(figName); 

    open(figName);
    hFig = gcf();
    axesObjs = get(hFig, 'Children');  %axes handles
    tmp = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
    dataObjs = tmp{2,1};

    res(:,1) = get(dataObjs(2), 'XData');  %data from low-level grahics objects
    res(:,2) = get(dataObjs(2), 'YData')-get(dataObjs(6), 'YData');
    res(:,3) = get(dataObjs(7), 'YData')-get(dataObjs(6), 'YData');
    
    norm = max(res(:,3));
    res(:,3) = res(:,3)/norm;
    res(:,2) = res(:,2)/norm;
    
    figure(2);
       plot(res(:,1),res(:,2),'-g',res(:,1),res(:,3),'or');
       
    csvwrite(strcat(name,'.csv'), res);
