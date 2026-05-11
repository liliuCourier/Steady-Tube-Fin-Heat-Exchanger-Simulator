%% tubeConnectionDesigner_final_robust.m
% 单节点显示，支持总入口和总出口，允许合流，支持单条连线删除
% 总出口平面由上游决定（不预设出口侧）
% 导出矩阵格式：行i，列j：-1上游，1下游
% 增强版：重置功能可重复使用，彻底清除连线图形
% 新增：进口方向设置，出口流向一致性检查
function  HX_Path_Planner1()

%% 获取用户输入
prompt = {'Tube per bank:', ...
          'Bank number:', ...
          'Inlet number:', ...
          'Outlet number:'};
dlgtitle = '管束参数输入';
dims = [1 40];
definput = {'3', '2', '1', '1'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    return;
end

rows = str2double(answer{1});
cols = str2double(answer{2});
numInlets = str2double(answer{3});
numOutlets = str2double(answer{4});

if any(isnan([rows, cols, numInlets, numOutlets])) || ...
   any([rows, cols, numInlets, numOutlets] < 1) || ...
   any(mod([rows, cols, numInlets, numOutlets],1) ~= 0)
    errordlg('所有输入必须为正整数', '输入错误');
    return;
end

%% 计算管节点坐标
nTubes = rows * cols;
x0 = 1:cols;
y0 = 1:rows;
[X, Y] = meshgrid(x0, y0);
tubeX = X(:)';
tubeY = Y(:)';

%% 创建图形界面
fig = figure('Name', 'Tube Connection Designer', 'NumberTitle', 'off', ...
             'Position', [100 100 1000 650], 'MenuBar', 'none');

% 左侧面板
leftPanel = uipanel('Parent', fig, 'Position', [0.02 0.15 0.18 0.75], 'Title', '状态');

% 入口状态
uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '总入口管:', ...
    "Units","normalized",'Position', [0 0.9 1 0.1], 'HorizontalAlignment', 'left');
txtInletList = uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '无', ...
    "Units","normalized",'Position', [0 0.8 1 0.1], 'HorizontalAlignment', 'left', 'ForegroundColor', 'r');

% 出口状态
uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '总出口管:', ...
    "Units","normalized",'Position', [0 0.7 1 0.1], 'HorizontalAlignment', 'left');
txtOutletList = uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '无', ...
    "Units","normalized",'Position', [0 0.6 1 0.1], 'HorizontalAlignment', 'left', 'ForegroundColor', 'b');

% 连线数量
uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '连线数量:', ...
    "Units","normalized",'Position', [0 0.5 1 0.1], 'HorizontalAlignment', 'left');
txtConnCount = uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '0', ...
    "Units","normalized",'Position', [0 0.4 1 0.1], 'HorizontalAlignment', 'left');

% 提示信息
txtPrompt = uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '', ...
    "Units","normalized",'Position', [0 0.2 1 0.1], 'HorizontalAlignment', 'left', 'ForegroundColor', 'k');

% 新增：进口方向设置
uicontrol('Parent', leftPanel, 'Style', 'text', 'String', '进口方向:', ...
    "Units","normalized",'Position', [0 0.05 1 0.1], 'HorizontalAlignment', 'left');
pmInletDir = uicontrol('Parent', leftPanel, 'Style', 'popupmenu', ...
    'String', {'左侧(前)', '右侧(后)'}, 'Value', 1, ...
    "Units","normalized",'Position', [0 0 0.5 0.1], 'Callback', @cbInletDir);

% 绘图区域
ax = axes('Parent', fig, 'Position', [0.25 0.2 0.7 0.7]);
 hold on; grid on;%axis equal;
xlabel('Bank'); ylabel('Bank Tube');
xlim([0.5, cols+0.5]); ylim([0.5, rows+0.5]);
xtickslabel = 1:cols;
ytickslabel = 1:rows;
xticks(xtickslabel)
yticks(ytickslabel)
set(ax, 'YDir', 'reverse');

% 绘制所有管节点
hTubes = gobjects(1, nTubes);
for i = 1:nTubes
    hTubes(i) = plot(tubeX(i), tubeY(i), 'o', ...
        'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5, 'MarkerSize', 12, 'Tag', 'tube');
end

% 显示管编号
for i = 1:nTubes
    text(tubeX(i)+0.1, tubeY(i)-0.1, num2str(i), 'FontSize', 10, 'Color', 'k', 'Tag', 'tubeLabel');
end

title(sprintf('%dBanks × %dTubes（共%d根管）',cols, rows, nTubes));

bottom_location.xmin = 0.27;
bottom_location.ymin = 0.01;

bottom_width = 0.15;
bottom_height = 0.05;

bottom_y_spacing  = 0.065;
bottom_x_spacing  = 0.17;


% 按钮
btnInlet = uicontrol('Style', 'pushbutton', 'String', sprintf('设置入口 (需%d个)', numInlets), ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin bottom_location.ymin + bottom_y_spacing  bottom_width bottom_height], 'Callback', @cbInlet);
btnOutlet = uicontrol('Style', 'pushbutton', 'String', sprintf('设置出口 (需%d个)', numOutlets), ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + bottom_x_spacing bottom_location.ymin + bottom_y_spacing bottom_width bottom_height], 'Callback', @cbOutlet);
btnConnect = uicontrol('Style', 'togglebutton', 'String', '连线模式', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + 2*bottom_x_spacing bottom_location.ymin + bottom_y_spacing bottom_width bottom_height], 'Callback', @cbConnect);
btnDelete = uicontrol('Style', 'togglebutton', 'String', '删除模式', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + 3*bottom_x_spacing bottom_location.ymin + bottom_y_spacing bottom_width bottom_height], 'Callback', @cbDelete);
btnResetConn = uicontrol('Style', 'pushbutton', 'String', '重置所有连线', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin bottom_location.ymin bottom_width bottom_height], 'Callback', @cbResetConn);
btnClear = uicontrol('Style', 'pushbutton', 'String', '清除进出口', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + 1*bottom_x_spacing bottom_location.ymin bottom_width bottom_height], 'Callback', @cbClear);
btnExport = uicontrol('Style', 'pushbutton', 'String', '导出结果', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + 2*bottom_x_spacing bottom_location.ymin bottom_width bottom_height], 'Callback', @cbExport);
btnRestart = uicontrol('Style', 'pushbutton', 'String', '重新输入', ...
    "Units","normalized","Parent",fig,...
    'Position', [bottom_location.xmin + 3*bottom_x_spacing bottom_location.ymin bottom_width bottom_height], 'Callback', @cbRestart);

% 底部状态条
txtStatus = uicontrol('Style', 'text', 'String', '就绪', ...
    "Units","normalized","Parent",fig,...
    'Position', [0.05 0 0.1 0.10], 'HorizontalAlignment', 'left');

%% 初始化数据
handles = struct();
handles.fig = fig;
handles.ax = ax;
handles.tubeX = tubeX;
handles.tubeY = tubeY;
handles.nTubes = nTubes;
handles.hTubes = hTubes;
handles.txtInletList = txtInletList;
handles.txtOutletList = txtOutletList;
handles.txtConnCount = txtConnCount;
handles.txtPrompt = txtPrompt;
handles.txtStatus = txtStatus;
handles.btnInlet = btnInlet;
handles.btnOutlet = btnOutlet;
handles.btnConnect = btnConnect;
handles.btnDelete = btnDelete;
handles.btnResetConn = btnResetConn;
handles.btnClear = btnClear;
handles.pmInletDir = pmInletDir;    % 新增

handles.numInlets = numInlets;
handles.numOutlets = numOutlets;
handles.selectedInlets = [];
handles.selectedOutlets = [];
handles.recordingMode = '';      % 'inlet', 'outlet', 'connect', 'delete'
handles.remaining = 0;
handles.inletDirection = 1;      % 新增：1=左进(前), 2=右进(后)

handles.tubeInletSide = zeros(1, nTubes);   % 1左 2右 0未定义
handles.tubeOutletSide = zeros(1, nTubes);  % 1左 2右 0未定义（总出口此值保持0）

handles.connections = struct('startTube', {}, 'endTube', {}, 'lineHandle', {});
handles.adjList = cell(1, nTubes);
handles.row = rows;
handles.col = cols;

handles.connectFirstTube = [];

guidata(fig, handles);

set(fig, 'WindowButtonDownFcn', @figClick);
set(fig, 'WindowButtonMotionFcn', @figMove);

updateDisplay(handles);


end

%% 局部函数定义
function updateDisplay(handles)
hTubes = handles.hTubes;
inlets = handles.selectedInlets;
outlets = handles.selectedOutlets;
txtInletList = handles.txtInletList;
txtOutletList = handles.txtOutletList;
txtConnCount = handles.txtConnCount;
txtPrompt = handles.txtPrompt;
recordingMode = handles.recordingMode;
remaining = handles.remaining;
connections = handles.connections;

for i = 1:length(hTubes)
    set(hTubes(i), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none', 'Marker', 'o');
end

for idx = inlets
    set(hTubes(idx), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'o');
end

for idx = outlets
    set(hTubes(idx), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Marker', 'o');
end

inletStrs = {};
for i = 1:length(inlets)
    t = inlets(i);
    side = handles.tubeInletSide(t);
    if side == 1
        inletStrs{end+1} = sprintf('%d(左)', t);
    elseif side == 2
        inletStrs{end+1} = sprintf('%d(右)', t);
    else
        inletStrs{end+1} = sprintf('%d(?)', t);
    end
end
set(txtInletList, 'String', iff(isempty(inletStrs), '无', strjoin(inletStrs, ', ')));

if isempty(outlets)
    set(txtOutletList, 'String', '无');
else
    outletStrs = arrayfun(@(x) num2str(x), outlets, 'UniformOutput', false);
    set(txtOutletList, 'String', strjoin(outletStrs, ', '));
end

set(txtConnCount, 'String', num2str(length(connections)));

if strcmp(recordingMode, 'inlet')
    set(txtPrompt, 'String', sprintf('正在选择总入口管，还需 %d 个', remaining));
elseif strcmp(recordingMode, 'outlet')
    set(txtPrompt, 'String', sprintf('正在选择总出口管，还需 %d 个', remaining));
elseif get(handles.btnConnect, 'Value') == 1
    if isempty(handles.connectFirstTube)
        set(txtPrompt, 'String', '连线模式：点击起点管');
    else
        set(txtPrompt, 'String', sprintf('连线模式：已选起点管 %d，点击终点管', handles.connectFirstTube));
    end
elseif get(handles.btnDelete, 'Value') == 1
    set(txtPrompt, 'String', '删除模式：点击连线即可删除');
else
    set(txtPrompt, 'String', '');
end

guidata(handles.fig, handles);
end

function cbInletDir(hObject, ~)
handles = guidata(hObject);
newDir = get(hObject, 'Value');
oldDir = handles.inletDirection;
if newDir == oldDir
    return;   % 方向未变，无需操作
end
handles.inletDirection = newDir;

% 如果已经存在任何进口管或内部连线，则全局翻转所有管的侧信息
if ~isempty(handles.selectedInlets) || ~isempty(handles.connections)
    % 辅助函数：1变2，2变1，0不变
    flipSide = @(s) iff(s==1, 2, iff(s==2, 1, 0));
    for i = 1:handles.nTubes
        handles.tubeInletSide(i)  = flipSide(handles.tubeInletSide(i));
        handles.tubeOutletSide(i) = flipSide(handles.tubeOutletSide(i));
    end
else
    % 没有任何连接，仅根据新方向调整已有进口管的侧信息
    for t = handles.selectedInlets
        if newDir == 1
            handles.tubeInletSide(t) = 1;
            handles.tubeOutletSide(t) = 2;
        else
            handles.tubeInletSide(t) = 2;
            handles.tubeOutletSide(t) = 1;
        end
    end
end

% 更新连线虚实和界面显示
drawAllConnections(handles);
updateDisplay(handles);
guidata(hObject, handles);
end

function drawAllConnections(handles)
delete(findobj(handles.ax, 'Tag', 'connection'));
newConns = handles.connections;
for k = 1:length(newConns)
    s = newConns(k).startTube;
    e = newConns(k).endTube;
    x1 = handles.tubeX(s); y1 = handles.tubeY(s);
    x2 = handles.tubeX(e); y2 = handles.tubeY(e);
    outSide = handles.tubeOutletSide(s);
    lineStyle = iff(outSide == 1, '-', '--');   % 1左(前)→实线，2右(后)→虚线
    h = quiver(handles.ax, x1, y1, x2-x1, y2-y1, 0, ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
        'LineStyle', lineStyle, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off', ...
        'Tag', 'connection', 'UserData', k);
    newConns(k).lineHandle = h;
end
handles.connections = newConns;
guidata(handles.fig, handles);
end

function recomputeSides(handles)
    nTubes = handles.nTubes;
    inletSide = zeros(1, nTubes);
    outletSide = zeros(1, nTubes);

    % 第一步：保留进口管的侧信息（根据当前进口方向）
    for t = handles.selectedInlets
        if handles.inletDirection == 1   % 左进
            inletSide(t) = 1;
            outletSide(t) = 2;
        else                            % 右进
            inletSide(t) = 2;
            outletSide(t) = 1;
        end
    end
    % 所有普通管和出口管的侧暂时为0（出口管也会由传播得到）

    % 第二步：通过连线传播侧信息
    changed = true;
    while changed
        changed = false;
        for k = 1:length(handles.connections)
            s = handles.connections(k).startTube;
            e = handles.connections(k).endTube;

            % 若起点已有入口侧，则可确定其出口侧
            if inletSide(s) ~= 0
                expectedOut = 3 - inletSide(s);
                if outletSide(s) == 0
                    outletSide(s) = expectedOut;
                    changed = true;
                elseif outletSide(s) ~= expectedOut
                    error('起点 %d 出口侧不一致', s);
                end
            end

            % 若起点已有出口侧，则可确定终点的入口侧
            if outletSide(s) ~= 0
                expectedIn = outletSide(s);
                if inletSide(e) == 0
                    inletSide(e) = expectedIn;
                    changed = true;
                elseif inletSide(e) ~= expectedIn
                    error('终点 %d 入口侧不一致', e);
                end
            end

            % 若终点已有入口侧，且不是出口管，则可确定其出口侧
            if inletSide(e) ~= 0 && ~ismember(e, handles.selectedOutlets)
                expectedOut = 3 - inletSide(e);
                if outletSide(e) == 0
                    outletSide(e) = expectedOut;
                    changed = true;
                elseif outletSide(e) ~= expectedOut
                    error('管 %d 出口侧不一致', e);
                end
            end
        end
    end

    handles.tubeInletSide = inletSide;
    handles.tubeOutletSide = outletSide;
    guidata(handles.fig, handles);
end

% function deleteConnection(handles, idx)
% if idx < 1 || idx > length(handles.connections)
%     return;
% end
% if ishandle(handles.connections(idx).lineHandle)
%     delete(handles.connections(idx).lineHandle);
% end
% s = handles.connections(idx).startTube;
% e = handles.connections(idx).endTube;
% handles.adjList{s}(handles.adjList{s} == e) = [];
% handles.connections(idx) = [];
% recomputeSides(handles);
% drawAllConnections(handles);
% updateDisplay(handles);
% end

function deleteConnection(handles, idx)
if idx < 1 || idx > length(handles.connections)
    return;
end
if ishandle(handles.connections(idx).lineHandle)
    delete(handles.connections(idx).lineHandle);
end
s = handles.connections(idx).startTube;
e = handles.connections(idx).endTube;
handles.adjList{s}(handles.adjList{s} == e) = [];
handles.connections(idx) = [];

% 重置受影响管的侧信息（仅限普通管，进出口管不动）
% 检查起点管是否还有其他出边
if ~ismember(s, handles.selectedOutlets)%~ismember(s, handles.selectedInlets) && 
    if ~any(arrayfun(@(c) c.startTube == s, handles.connections))
        handles.tubeOutletSide(s) = 0;   % 无出边，出口侧清零
    end
end
% 检查终点管是否还有其他入边
if ~ismember(e, handles.selectedInlets) %&& ~ismember(e, handles.selectedOutlets)
    if ~any(arrayfun(@(c) c.endTube == e, handles.connections))
        handles.tubeInletSide(e) = 0;   % 无入边，入口侧清零
    end
end

recomputeSides(handles);
drawAllConnections(handles);
updateDisplay(handles);
end

function cbInlet(hObject, ~)
handles = guidata(hObject);
if strcmp(handles.recordingMode, 'inlet')
    handles.recordingMode = '';
    handles.remaining = 0;
    set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets));
    set(handles.btnOutlet, 'Enable', 'on');
    set(handles.btnConnect, 'Enable', 'on');
    set(handles.btnDelete, 'Enable', 'on');
    set(handles.btnResetConn, 'Enable', 'on');
    set(handles.btnClear, 'Enable', 'on');
    set(handles.txtStatus, 'String', '已退出入口选择模式');
else
    if length(handles.selectedInlets) >= handles.numInlets
        set(handles.txtStatus, 'String', '入口已选满，如需修改请先清除');
        return;
    end
    if get(handles.btnConnect, 'Value') == 1
        set(handles.btnConnect, 'Value', 0);
        handles.connectFirstTube = [];
    end
    if get(handles.btnDelete, 'Value') == 1
        set(handles.btnDelete, 'Value', 0);
    end
    if strcmp(handles.recordingMode, 'outlet')
        handles.recordingMode = '';
        set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets));
    end
    handles.recordingMode = 'inlet';
    handles.remaining = handles.numInlets - length(handles.selectedInlets);
    set(handles.btnInlet, 'String', '取消');
    set(handles.btnOutlet, 'Enable', 'off');
    set(handles.btnConnect, 'Enable', 'off');
    set(handles.btnDelete, 'Enable', 'off');
    set(handles.btnResetConn, 'Enable', 'off');
    set(handles.btnClear, 'Enable', 'off');
    set(handles.txtStatus, 'String', sprintf('请点击 %d 个管作为总入口', handles.remaining));
end
guidata(hObject, handles);
updateDisplay(handles);
end

function cbOutlet(hObject, ~)
handles = guidata(hObject);
if strcmp(handles.recordingMode, 'outlet')
    handles.recordingMode = '';
    handles.remaining = 0;
    set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets));
    set(handles.btnInlet, 'Enable', 'on');
    set(handles.btnConnect, 'Enable', 'on');
    set(handles.btnDelete, 'Enable', 'on');
    set(handles.btnResetConn, 'Enable', 'on');
    set(handles.btnClear, 'Enable', 'on');
    set(handles.txtStatus, 'String', '已退出出口选择模式');
else
    if length(handles.selectedOutlets) >= handles.numOutlets
        set(handles.txtStatus, 'String', '出口已选满，如需修改请先清除');
        return;
    end
    if get(handles.btnConnect, 'Value') == 1
        set(handles.btnConnect, 'Value', 0);
        handles.connectFirstTube = [];
    end
    if get(handles.btnDelete, 'Value') == 1
        set(handles.btnDelete, 'Value', 0);
    end
    if strcmp(handles.recordingMode, 'inlet')
        handles.recordingMode = '';
        set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets));
    end
    handles.recordingMode = 'outlet';
    handles.remaining = handles.numOutlets - length(handles.selectedOutlets);
    set(handles.btnOutlet, 'String', '取消');
    set(handles.btnInlet, 'Enable', 'off');
    set(handles.btnConnect, 'Enable', 'off');
    set(handles.btnDelete, 'Enable', 'off');
    set(handles.btnResetConn, 'Enable', 'off');
    set(handles.btnClear, 'Enable', 'off');
    set(handles.txtStatus, 'String', sprintf('请点击 %d 个管作为总出口', handles.remaining));
end
guidata(hObject, handles);
updateDisplay(handles);
end

function cbConnect(hObject, ~)
handles = guidata(hObject);
if get(hObject, 'Value') == 1
    if ~isempty(handles.recordingMode)
        handles.recordingMode = '';
        handles.remaining = 0;
        set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets));
        set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets));
    end
    if get(handles.btnDelete, 'Value') == 1
        set(handles.btnDelete, 'Value', 0);
    end
    handles.connectFirstTube = [];
    set(handles.btnInlet, 'Enable', 'off');
    set(handles.btnOutlet, 'Enable', 'off');
    set(handles.btnDelete, 'Enable', 'off');
    set(handles.btnResetConn, 'Enable', 'off');
    set(handles.btnClear, 'Enable', 'off');
    set(handles.txtStatus, 'String', '连线模式已开启，点击起点管');
else
    handles.connectFirstTube = [];
    set(handles.btnInlet, 'Enable', 'on');
    set(handles.btnOutlet, 'Enable', 'on');
    set(handles.btnDelete, 'Enable', 'on');
    set(handles.btnResetConn, 'Enable', 'on');
    set(handles.btnClear, 'Enable', 'on');
    set(handles.txtStatus, 'String', '已退出连线模式');
end
guidata(hObject, handles);
updateDisplay(handles);
end

function cbDelete(hObject, ~)
handles = guidata(hObject);
if get(hObject, 'Value') == 1
    if ~isempty(handles.recordingMode)
        handles.recordingMode = '';
        handles.remaining = 0;
        set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets));
        set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets));
    end
    if get(handles.btnConnect, 'Value') == 1
        set(handles.btnConnect, 'Value', 0);
        handles.connectFirstTube = [];
    end
    set(handles.btnInlet, 'Enable', 'off');
    set(handles.btnOutlet, 'Enable', 'off');
    set(handles.btnConnect, 'Enable', 'off');
    set(handles.btnResetConn, 'Enable', 'off');
    set(handles.btnClear, 'Enable', 'off');
    set(handles.txtStatus, 'String', '删除模式已开启，点击连线删除');
else
    set(handles.btnInlet, 'Enable', 'on');
    set(handles.btnOutlet, 'Enable', 'on');
    set(handles.btnConnect, 'Enable', 'on');
    set(handles.btnResetConn, 'Enable', 'on');
    set(handles.btnClear, 'Enable', 'on');
    set(handles.txtStatus, 'String', '已退出删除模式');
end
guidata(hObject, handles);
updateDisplay(handles);
end

function cbResetConn(hObject, ~)
handles = guidata(hObject);
delete(findobj(handles.ax, 'Tag', 'connection'));
handles.connections = struct('startTube', {}, 'endTube', {}, 'lineHandle', {});
handles.adjList = cell(1, handles.nTubes);
% 重置侧信息：出口管侧清空，进口管根据当前方向重新设定
for t = 1:handles.nTubes
    if ismember(t, handles.selectedInlets)
        if handles.inletDirection == 1
            handles.tubeInletSide(t) = 1;
            handles.tubeOutletSide(t) = 2;
        else
            handles.tubeInletSide(t) = 2;
            handles.tubeOutletSide(t) = 1;
        end
    elseif ismember(tcbIn, handles.selectedOutlets)
        handles.tubeInletSide(t) = 0;
        handles.tubeOutletSide(t) = 0;
    else
        handles.tubeInletSide(t) = 0;
        handles.tubeOutletSide(t) = 0;
    end
end
handles.connectFirstTube = [];
set(handles.txtStatus, 'String', '所有连线已清除');
guidata(hObject, handles);
updateDisplay(handles);
end

function cbClear(hObject, ~)
handles = guidata(hObject);
handles.selectedInlets = [];
handles.selectedOutlets = [];
handles.recordingMode = '';
handles.remaining = 0;
handles.tubeInletSide = zeros(1, handles.nTubes);
handles.tubeOutletSide = zeros(1, handles.nTubes);
delete(findobj(handles.ax, 'Tag', 'connection'));
handles.connections = struct('startTube', {}, 'endTube', {}, 'lineHandle', {});
handles.adjList = cell(1, handles.nTubes);
handles.connectFirstTube = [];
set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets), 'Enable', 'on');
set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets), 'Enable', 'on');
set(handles.btnConnect, 'Enable', 'on');
set(handles.btnDelete, 'Enable', 'on');
set(handles.btnResetConn, 'Enable', 'on');
set(handles.btnClear, 'Enable', 'on');
set(handles.txtStatus, 'String', '已清除所有进出口和连线');
guidata(hObject, handles);
updateDisplay(handles);
end

function cbExport(hObject, ~)
handles = guidata(hObject);
nTubes = handles.nTubes;
edges = handles.connections;
inlets = handles.selectedInlets;
outlets = handles.selectedOutlets;
row = handles.row;
col = handles.col;

% 确保侧信息是最新的
recomputeSides(handles);
handles = guidata(hObject);  % 重新获取更新后的 handles

% ==================== 完整性检查 ====================
outDegree = zeros(1, nTubes);
inDegree = zeros(1, nTubes);
for k = 1:length(edges)
    s = edges(k).startTube;
    e = edges(k).endTube;
    outDegree(s) = outDegree(s) + 1;
    inDegree(e) = inDegree(e) + 1;
end

errMsg = '';
for t = inlets
    if outDegree(t) == 0
        errMsg = [errMsg, sprintf('进口管 %d 没有下游连接！\n', t)];
    end
    if inDegree(t) > 0
        errMsg = [errMsg, sprintf('进口管 %d 不应有上游连接！\n', t)];
    end
end
for t = outlets
    if inDegree(t) == 0
        errMsg = [errMsg, sprintf('出口管 %d 没有上游连接！\n', t)];
    end
    if outDegree(t) > 0
        errMsg = [errMsg, sprintf('出口管 %d 不应有下游连接！\n', t)];
    end
end
for t = 1:nTubes
    if ismember(t, inlets) || ismember(t, outlets)
        continue;
    end
    if inDegree(t) == 0
        errMsg = [errMsg, sprintf('管 %d 缺少上游连接（流体无法流入）！\n', t)];
    end
    if outDegree(t) == 0
        errMsg = [errMsg, sprintf('管 %d 缺少下游连接（流体无法流出）！\n', t)];
    end  
end

% 可达性检查
adj = cell(1, nTubes);
for k = 1:length(edges)
    s = edges(k).startTube;
    e = edges(k).endTube;
    adj{s} = [adj{s}, e];
end
visited = false(1, nTubes);
stack = inlets(:)';
for t = inlets
    visited(t) = true;
end
while ~isempty(stack)
    u = stack(end); stack(end) = [];
    for v = adj{u}
        if ~visited(v)
            visited(v) = true;
            stack = [stack, v];
        end
    end
end
unreachable = find(~visited);
if ~isempty(unreachable)
    errMsg = [errMsg, sprintf('以下管无法从进口到达：%s\n', num2str(unreachable))];
end
for t = outlets
    if ~visited(t)
        errMsg = [errMsg, sprintf('出口管 %d 无法从进口到达！\n', t)];
    end
end

if ~isempty(errMsg)
    errordlg(errMsg, '流路不完整，无法导出');
    set(handles.txtStatus, 'String', '导出失败：流路不完整，请检查连接');
    return;
end

% ==================== 出口流向一致性警告 ====================
% ==================== 构建关联矩阵 ====================
nodeMap = containers.Map();
nextNode = 1;

% 辅助函数（嵌套）
    function nid = getOrCreateNode(tube, side)
        key = sprintf('%d_%s', tube, side);
        if isKey(nodeMap, key)
            nid = nodeMap(key);
        else
            nid = nextNode;
            nodeMap(key) = nid;
            nextNode = nextNode + 1;
        end
    end

% 进口联箱
inletHeaderNode = nextNode; nextNode = nextNode + 1;
for t = inlets
    key_in = sprintf('%d_in', t);
    nodeMap(key_in) = inletHeaderNode;
end
% 出口联箱
outletHeaderNode = nextNode; nextNode = nextNode + 1;
for t = outlets
    key_out = sprintf('%d_out', t);
    nodeMap(key_out) = outletHeaderNode;
end

% 内部U型弯连接
for k = 1:length(edges)
    s = edges(k).startTube;
    e = edges(k).endTube;
    key_so = sprintf('%d_out', s);
    key_ei = sprintf('%d_in', e);
    if isKey(nodeMap, key_so) && isKey(nodeMap, key_ei)
        nid_so = nodeMap(key_so);
        nid_ei = nodeMap(key_ei);
        if nid_so ~= nid_ei
            allKeys = keys(nodeMap);
            for i = 1:length(allKeys)
                if nodeMap(allKeys{i}) == nid_ei
                    nodeMap(allKeys{i}) = nid_so;
                end
            end
        end
    elseif isKey(nodeMap, key_so)
        nodeMap(key_ei) = nodeMap(key_so);
    elseif isKey(nodeMap, key_ei)
        nodeMap(key_so) = nodeMap(key_ei);
    else
        nid = nextNode; nextNode = nextNode + 1;
        nodeMap(key_so) = nid;
        nodeMap(key_ei) = nid;
    end
end

% 补全未分配端点
for i = 1:nTubes
    key_in = sprintf('%d_in', i);
    key_out = sprintf('%d_out', i);
    if ~isKey(nodeMap, key_in)
        nodeMap(key_in) = nextNode; nextNode = nextNode + 1;
    end
    if ~isKey(nodeMap, key_out)
        nodeMap(key_out) = nextNode; nextNode = nextNode + 1;
    end
end

numNodes = nextNode - 1;
M_full = zeros(numNodes, nTubes);
for tube = 1:nTubes
    key_in  = sprintf('%d_in',  tube);
    key_out = sprintf('%d_out', tube);
    n_in  = nodeMap(key_in);
    n_out = nodeMap(key_out);
    M_full(n_in,  tube) = -1;
    M_full(n_out, tube) =  1;
end

nodeLabels = cell(numNodes, 1);
for i = 1:numNodes
    nodeLabels{i} = sprintf('N%d', i);
end
nodeLabels{inletHeaderNode}  = 'INLET_HEADER';
nodeLabels{outletHeaderNode} = 'OUTLET_HEADER';

M_internal = M_full;
M_internal([inletHeaderNode, outletHeaderNode], :) = [];

% 流向信息
tubeFlowDir = zeros(1, nTubes);
for t = 1:nTubes
    if handles.tubeInletSide(t) == 1
        tubeFlowDir(t) = 1;
    elseif handles.tubeInletSide(t) == 2
        tubeFlowDir(t) = 0;
    else
        tubeFlowDir(t) = NaN;
    end
end
inletSide  = handles.tubeInletSide;
outletSide = handles.tubeOutletSide;

outflowdir = tubeFlowDir(outlets);
if all(outflowdir == outflowdir(1))
    if outflowdir(1) == tubeFlowDir(inlets(1))
        warning('进出口管位置不在一边');
    end
else
    inconsistent = outlets(outflowdir ~= outflowdir(1));
    consistent = outlets(outflowdir == outflowdir(1));
    errMsg = [errMsg, '出口管位置不一致'];
    errordlg(errMsg);
    error(errMsg)
    % warning('以下管与第一个管方向不一致：');
    % disp(inconsistent);
    % warning('以下管与第一个管方向一致：');
    % disp(consistent);
end

% 预处理矩阵信息

TC_matrix = -1*M_internal;
FlowDirectionInf = find(tubeFlowDir == 0);
IO_inlet = zeros(size(TC_matrix));
IO_outlet = zeros(size(TC_matrix));
IO_inlet(TC_matrix == 1) = 1;
IO_outlet(TC_matrix == -1) = -1;

% 生成U型弯长度示意向量


% 打包导出
TCinf.TC_matrix        = TC_matrix;
TCinf.inlet_num        = inlets;
TCinf.outlet_num       = outlets;
TCinf.FlowDirection    = tubeFlowDir;
TCinf.FlowDirectionInf = FlowDirectionInf;
TCinf.row              = row;
TCinf.col              = col;
TCinf.con_num          = size(M_internal,1);
TCinf.Tube_num         = row*col;
TCinf.IO_inlet         = IO_inlet;
TCinf.IO_outlet        = IO_outlet;
TCinf.Uband_L          = Uband_L;

assignin('base', 'TCinf', TCinf);

disp('=== 导出结果 ===');
disp('节点-管道关联矩阵(含虚拟节点):'); disp(M_internal);
disp('管流向 (0=前→后, 1=后→前):'); disp(tubeFlowDir);

set(handles.txtStatus, 'String', '导出成功：矩阵与流向信息已保存至 TCinf');
guidata(hObject, handles);
end

function cbRestart(hObject, ~)
handles = guidata(hObject);
close(handles.fig);
evalin('base', 'tubeConnectionDesigner_final_robust');
end

function figClick(hObject, ~)
handles = guidata(hObject);
cp = get(handles.ax, 'CurrentPoint');
x = cp(1,1); y = cp(1,2);

if get(handles.btnDelete, 'Value') == 1
    threshold = 0.3;
    minDist = inf;
    closestIdx = [];
    for k = 1:length(handles.connections)
        s = handles.connections(k).startTube;
        e = handles.connections(k).endTube;
        x1 = handles.tubeX(s); y1 = handles.tubeY(s);
        x2 = handles.tubeX(e); y2 = handles.tubeY(e);
        d = pointToLineDistance(x, y, x1, y1, x2, y2);
        if d < threshold && d < minDist
            minDist = d;
            closestIdx = k;
        end
    end
    if ~isempty(closestIdx)
        deleteConnection(handles, closestIdx);
    end
    return;
end

dist = sqrt((handles.tubeX - x).^2 + (handles.tubeY - y).^2);
[minDist, idx] = min(dist);
if minDist > 0.3
    return;
end
tube = idx;

if strcmp(handles.recordingMode, 'inlet')
    if ismember(tube, handles.selectedInlets)
        set(handles.txtStatus, 'String', sprintf('管 %d 已是总入口', tube));
        return;
    end
    if ismember(tube, handles.selectedOutlets)
        set(handles.txtStatus, 'String', sprintf('管 %d 已是总出口，不能作为总入口', tube));
        return;
    end
    for k = 1:length(handles.connections)
        if handles.connections(k).endTube == tube
            set(handles.txtStatus, 'String', sprintf('管 %d 已有上游连线，不能设为总入口', tube));
            return;
        end
    end
    handles.selectedInlets(end+1) = tube;
    % 设置进口管侧信息（根据当前进口方向）
    if handles.inletDirection == 1
        handles.tubeInletSide(tube) = 1;
        handles.tubeOutletSide(tube) = 2;
    else
        handles.tubeInletSide(tube) = 2;
        handles.tubeOutletSide(tube) = 1;
    end
    handles.remaining = handles.remaining - 1;
    set(handles.txtStatus, 'String', sprintf('已选总入口管 %d，还需 %d 个', tube, handles.remaining));
    if handles.remaining == 0
        handles.recordingMode = '';
        set(handles.btnInlet, 'String', sprintf('设置入口 (需%d个)', handles.numInlets), 'Enable', 'on');
        set(handles.btnOutlet, 'Enable', 'on');
        set(handles.btnConnect, 'Enable', 'on');
        set(handles.btnDelete, 'Enable', 'on');
        set(handles.btnResetConn, 'Enable', 'on');
        set(handles.btnClear, 'Enable', 'on');
        set(handles.txtStatus, 'String', '入口选择完成');
    end

elseif strcmp(handles.recordingMode, 'outlet')
    if ismember(tube, handles.selectedOutlets)
        set(handles.txtStatus, 'String', sprintf('管 %d 已是总出口', tube));
        return;
    end
    if ismember(tube, handles.selectedInlets)
        set(handles.txtStatus, 'String', sprintf('管 %d 已是总入口，不能作为总出口', tube));
        return;
    end
    for k = 1:length(handles.connections)
        if handles.connections(k).startTube == tube
            set(handles.txtStatus, 'String', sprintf('管 %d 已有下游连线，不能设为总出口', tube));
            return;
        end
    end
    handles.selectedOutlets(end+1) = tube;
    % 出口管的侧信息暂不设置，将由连接决定
    handles.tubeOutletSide(tube) = 0;
    handles.tubeInletSide(tube) = 0;
    handles.remaining = handles.remaining - 1;
    set(handles.txtStatus, 'String', sprintf('已选总出口管 %d，还需 %d 个', tube, handles.remaining));
    if handles.remaining == 0
        handles.recordingMode = '';
        set(handles.btnOutlet, 'String', sprintf('设置出口 (需%d个)', handles.numOutlets), 'Enable', 'on');
        set(handles.btnInlet, 'Enable', 'on');
        set(handles.btnConnect, 'Enable', 'on');
        set(handles.btnDelete, 'Enable', 'on');
        set(handles.btnResetConn, 'Enable', 'on');
        set(handles.btnClear, 'Enable', 'on');
        set(handles.txtStatus, 'String', '出口选择完成');
    end

elseif get(handles.btnConnect, 'Value') == 1
    if isempty(handles.connectFirstTube)
        if handles.tubeInletSide(tube) == 0
            set(handles.txtStatus, 'String', sprintf('管 %d 尚未有入口侧，不能作为起点', tube));
            return;
        end
        if ismember(tube, handles.selectedOutlets)
            set(handles.txtStatus, 'String', sprintf('管 %d 是总出口，不能作为起点', tube));
            return;
        end
        handles.connectFirstTube = tube;
        set(handles.txtStatus, 'String', sprintf('已选起点管 %d，点击终点管', tube));
    else
        startTube = handles.connectFirstTube;
        endTube = tube;
        if startTube == endTube
            set(handles.txtStatus, 'String', '起点和终点不能相同，请重新选择起点');
            handles.connectFirstTube = [];
            return;
        end

        startOutSide = 3 - handles.tubeInletSide(startTube);
        endInSide = startOutSide;

        if ismember(endTube, handles.selectedInlets)
            set(handles.txtStatus, 'String', sprintf('管 %d 是总入口，不能作为终点', endTube));
            handles.connectFirstTube = [];
            return;
        end

        if handles.tubeInletSide(endTube) ~= 0
            if handles.tubeInletSide(endTube) ~= endInSide
                set(handles.txtStatus, 'String', sprintf('管 %d 已有入口侧 %d，与当前所需 %d 不匹配', ...
                    endTube, handles.tubeInletSide(endTube), endInSide));
                handles.connectFirstTube = [];
                return;
            end
        end

        if ismember(endTube, handles.adjList{startTube})
            set(handles.txtStatus, 'String', '该连线已存在，请重新选择起点');
            handles.connectFirstTube = [];
            return;
        end

        tempAdj = handles.adjList;
        tempAdj{startTube} = [tempAdj{startTube}, endTube];
        if hasCycle(tempAdj, handles.nTubes)
            set(handles.txtStatus, 'String', '错误：该连线会导致环路，已取消');
            handles.connectFirstTube = [];
            return;
        end

        newConn = struct('startTube', startTube, 'endTube', endTube, 'lineHandle', []);
        x1 = handles.tubeX(startTube); y1 = handles.tubeY(startTube);
        x2 = handles.tubeX(endTube); y2 = handles.tubeY(endTube);
        lineStyle = iff(startOutSide == 1, '-', '--');
        h = quiver(handles.ax, x1, y1, x2-x1, y2-y1, 0, ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
            'LineStyle', lineStyle, ...
            'MaxHeadSize', 0.3, 'AutoScale', 'off', ...
            'Tag', 'connection', 'UserData', length(handles.connections)+1);
        newConn.lineHandle = h;

        handles.connections = [handles.connections, newConn];
        handles.adjList{startTube} = [handles.adjList{startTube}, endTube];

        handles.tubeOutletSide(startTube) = startOutSide;
        if handles.tubeInletSide(endTube) == 0
            handles.tubeInletSide(endTube) = endInSide;
        end
        if ~ismember(endTube, handles.selectedOutlets) && handles.tubeOutletSide(endTube) == 0
            handles.tubeOutletSide(endTube) = 3 - endInSide;
        end

        set(handles.txtStatus, 'String', sprintf('已创建连线 %d -> %d (%s)', ...
            startTube, endTube, iff(lineStyle=='-','实线','虚线')));
        handles.connectFirstTube = [];
    end
else
    set(handles.txtStatus, 'String', sprintf('点击管 %d', tube));
    return;
end

guidata(hObject, handles);
updateDisplay(handles);
drawAllConnections(handles);
end

function figMove(src, ~)
handles = guidata(src);
if isempty(handles) || isempty(handles.connections)
    return;
end

cp = get(handles.ax, 'CurrentPoint');
x = cp(1,1); y = cp(1,2);
threshold = 0.2;
minDist = inf;
closestIdx = [];

for k = 1:length(handles.connections)
    s = handles.connections(k).startTube;
    e = handles.connections(k).endTube;
    x1 = handles.tubeX(s); y1 = handles.tubeY(s);
    x2 = handles.tubeX(e); y2 = handles.tubeY(e);
    d = pointToLineDistance(x, y, x1, y1, x2, y2);
    if d < threshold && d < minDist
        minDist = d;
        closestIdx = k;
    end
end

for k = 1:length(handles.connections)
    if ishandle(handles.connections(k).lineHandle)
        set(handles.connections(k).lineHandle, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    end
end

if ~isempty(closestIdx)
    h = handles.connections(closestIdx).lineHandle;
    if ishandle(h)
        set(h, 'Color', 'k', 'LineWidth', 2.5);
        set(handles.txtStatus, 'String', sprintf('连线 %d -> %d', ...
            handles.connections(closestIdx).startTube, handles.connections(closestIdx).endTube));
    end
end
end

function d = pointToLineDistance(px, py, x1, y1, x2, y2)
vx = x2 - x1; vy = y2 - y1;
wx = px - x1; wy = py - y1;
c1 = wx*vx + wy*vy;
if c1 <= 0
    d = sqrt((px-x1)^2 + (py-y1)^2);
    return;
end
c2 = vx*vx + vy*vy;
if c2 <= c1
    d = sqrt((px-x2)^2 + (py-y2)^2);
    return;
end
b = c1 / c2;
pbx = x1 + b*vx;
pby = y1 + b*vy;
d = sqrt((px-pbx)^2 + (py-pby)^2);
end

function cycle = hasCycle(adjList, nNodes)
color = zeros(1, nNodes);
cycle = false;
for v = 1:nNodes
    if color(v) == 0
        if dfs(v)
            cycle = true;
            return;
        end
    end
end

    function found = dfs(u)
        color(u) = 1;
        for w = adjList{u}
            if color(w) == 1
                found = true;
                return;
            elseif color(w) == 0
                if dfs(w)
                    found = true;
                    return;
                end
            end
        end
        color(u) = 2;
        found = false;
    end
end

function s = iff(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end


