% PreProcessing — 预处理集成入口 (Demo1.02)
% 运行此脚本完成流路设计、几何设置、边界条件设置
% 完成后运行 Main 即可开始求解

% 添加子模块路径
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'PreProc'));
addpath(fullfile(root, 'Lib'));
addpath(fullfile(root, 'Solver'));

%% 0. 物性加载（如未加载）
if ~exist('Prop_handle', 'var')
    refprop_location = 'E:\refprop10\REFPROP';
    R = 'R134a';
    Prop_handle = Prop_load(refprop_location, R, 1e-3, 5.5, 80, 510, 100, 25, 25);
    fprintf('Prop_handle 已加载 (工质: %s)\n', R);
end

%% 1. 流路设计
fprintf('\n========== 步骤 1/3: 流路设计 ==========\n');
fprintf('请在打开的流路设计窗口中完成以下操作:\n');
fprintf('  1) 设置进出口管\n');
fprintf('  2) 绘制 U 型弯连线\n');
fprintf('  3) 点击 [导出结果] 按钮\n');
fprintf('设计窗口关闭按钮将隐藏窗口而非销毁，可随时通过 figure(hxDesigner) 唤出。\n');

HX_Path_Planner1();

fprintf('\n完成流路设计并导出后，按任意键继续...\n');
pause;

if ~exist('TCinf', 'var')
    error('未检测到 TCinf，请先在流路设计窗口中点击 [导出结果]');
end
fprintf('TCinf 已就绪 (%d管, %d排×%d列)\n', TCinf.Tube_num, TCinf.row, TCinf.col);

%% 2. 几何参数设置
fprintf('\n========== 步骤 2/3: 几何参数设置 ==========\n');
GeoCondition = GenerateGeo(TCinf);
assignin('base', 'GeoCondition', GeoCondition);
fprintf('GeoCondition 已保存至工作区\n');

%% 3. 边界条件设置
fprintf('\n========== 步骤 3/3: 边界条件设置 ==========\n');
BDCondition = GenerateBD(GeoCondition, Prop_handle);
assignin('base', 'BDCondition', BDCondition);
fprintf('BDCondition 已保存至工作区\n');

fprintf('\n========== 预处理完成 ==========\n');
fprintf('TCinf, GeoCondition, BDCondition 已就绪。\n');
fprintf('运行 Main 开始求解，求解完成后自动执行后处理。\n');
