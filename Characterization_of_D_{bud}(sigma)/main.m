A11=[-2,1;1,-3];
b11_temp=5;
A12=[-2,-1;-1,-10];
b12_temp=15;
A21=[-4,1;1,-4];
b21_temp=0;
A22=[-5,-1;-1,-2];
b22_temp=0;
x11_temp=3;
x12_temp=12;
x21_temp=25;
x22_temp=17;


[temp1,temp2]=testfun1([A11(1,:),b11_temp],A11(2,:),x11_temp);
x11=[x11_temp,temp1];
b11=[b11_temp;temp2];
[temp1,temp2]=testfun1([A12(1,:),b12_temp],A12(2,:),x12_temp);
x12=[x12_temp,temp1];
b12=[b12_temp;temp2];
[temp1,temp2]=testfun1([A21(2,:),b21_temp],A21(1,:),x21_temp);
x21=[x21_temp,temp1];
b21=[temp2;b21_temp];
[temp1,temp2]=testfun1([A22(2,:),b22_temp],A22(1,:),x22_temp);
x22=[x22_temp,temp1];
b22=[temp2;b22_temp];
%max point for sum J11+J21


% eta11=1;
% eta12=3;
% eta21=1;
% eta22=3;

eta11=1;
eta12=0;
eta21=5;
eta22=0;

lambda1=1;
lambda2=1;

UA=eta11*A11+eta21*A21+eta12*A12+eta22*A22;
Ub=eta11*b11+eta21*b21+eta12*b12+eta22*b22;

JA=A11+lambda1*A12+A21+lambda2*A22;
Jb=b11+lambda1*b12+b21+lambda2*b22;

xp1=-inv(JA)*(Jb);

xpe1=-inv(UA)*(Ub);
sx=5;
sy=5;


hFig = figure('Name', 'Interactive Sigma Adjustment', 'NumberTitle', 'off');

sigma=.5;


hScatter1 =scatter([xpe1(1)],[xpe1(2)],'filled', 'DisplayName', 'Jsum maximum Point');
hold on
hScatter2 = scatter([sx], [sy], 'filled', 'DisplayName', 'x0');
hImplicit = fimplicit(@(x,y) ([x,y]*JA*[x;y]/2+Jb'*[x;y])-([sx,sy]*JA*[sx;sy]/2+Jb'*[sx;sy]),[-50 50 -50 50],'r','DisplayName', 'Boundary of $\mathcal{D}(J_{\rm sum},x_0)$');
%fimplicit(@(x,y) ([x,y]*A11*[x;y]/2+b11'*[x;y]+c1)+([x,y]*A21*[x;y]/2+b21'*[x;y]+c2)-(xpe1'*A11*xpe1/2+b11'*xpe1+c1)-(xpe1'*A21*xpe1/2+b21'*xpe1+c2),[-50 50 -50 50],'r--');

fimplicit(@(x,y) ([x,y]*UA*[x;y]/2+Ub'*[x;y])-([sx,sy]*UA*[sx;sy]/2+Ub'*[sx;sy]),[-50 50 -50 50],'b', 'DisplayName', 'Boundary of $\mathcal{D}(\tilde{J}_{\rm sum},x_0)$');
Dbud=fimplicit(@(x,y) sigma*(([x,y]*UA*[x;y]/2+Ub'*[x;y])-([sx,sy]*UA*[sx;sy]/2+Ub'*[sx;sy]))-(([x,y]*JA*[x;y]/2+Jb'*[x;y])-([sx,sy]*JA*[sx;sy]/2+Jb'*[sx;sy])),[-100 100 -100 100],'m--', 'DisplayName', 'Boundary of $\mathcal{D}_{\rm bud}(\sigma)$');

% 将句柄存储在 guidata 中，以便在回调函数中更新

legend('show');  % 启动图例
legend('Location', 'best');  % 自动放置在最合适的位置

%guidata(hFig, struct('Dbud', Dbud));


hText = uicontrol('Style', 'text', 'String', 'Sigma: 0.5', ...  % 初始值Sigma为0.5
    'Units', 'normalized', 'Position', [0.1, 0.85, 0.8, 0.05], ...
    'FontSize', 12);

% --- 创建滑块控件后，将 hText 和其他对象一同存入 guidata ---
guidata(hFig, struct('Dbud', Dbud, 'hText', hText));

hSlider = uicontrol('Style', 'slider', ...
    'Min', 0.3, 'Max', 0.6, 'Value', 0.5, ...  % 初始值 0.5
    'SliderStep', [0.005/0.3, 0.1/0.3], ...
    'Units', 'normalized', 'Position', [0.1 0.02 0.8 0.05], ... % 滑条大小和位置
    'Callback', @(src, event) sliderCallback(src, hFig,hText, UA,Ub,JA,Jb, sx, sy)); 



% 定义滑条的回调函数
function sliderCallback(src, hFig,hText, UA,Ub,JA,Jb, sx, sy)
    % 获取当前的 sigma 值
    sigma = get(src, 'Value');

    data = guidata(hFig);
    if isfield(data, 'Dbud') && isvalid(data.Dbud)
        delete(data.Dbud);  % 删除旧 Dbud 曲线
    end
    hold on;
    % 根据新的 sigma，重新绘制 Dbud
    Dbud=fimplicit(@(x,y) sigma*(([x,y]*UA*[x;y]/2+Ub'*[x;y])-([sx,sy]*UA*[sx;sy]/2+Ub'*[sx;sy]))-(([x,y]*JA*[x;y]/2+Jb'*[x;y])-([sx,sy]*JA*[sx;sy]/2+Jb'*[sx;sy])),[-100 100 -100 100],'m--', 'DisplayName', 'Boundary of $\mathcal{D}_{\rm bud}(\sigma)$');
    % 保存新的 Dbud 句柄到 guidata 中，确保下次可以更新它
    data.Dbud = Dbud;
    guidata(hFig, data);

   if isfield(data, 'hText') && isvalid(data.hText)
        set(data.hText, 'String', sprintf('Sigma: %.3f', sigma));  % 更新显示
    else
        disp('hText 无效或未找到');
    end
end
