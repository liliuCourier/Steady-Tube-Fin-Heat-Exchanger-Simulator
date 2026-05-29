# Program_1_AI — 管翅式换热器稳态求解器

## 工作流
1. **PreProcessing.m**（预处理）— 含 GUI 交互，需用户手动操作：
   - HX_Path_Planner1（流路设计器，GUI）
   - GenerateGeo（几何参数对话框）
   - GenerateBD（边界条件对话框）
   - 完成后在 base workspace 生成 TCinf、GeoCondition、BDCondition、Prop_handle
2. **Main.m**（主求解器）— 依赖预处理变量，自动调用后处理
3. **PostProcessing.m**（后处理）— Main 完成后自动触发

## 重要约束
- 预处理（PreProcessing.m）包含用户交互界面，无法通过代码自动完成
- 每次运行 Main.m 前，必须先检查 TCinf、GeoCondition、BDCondition 是否已在 workspace 中存在
- 所有操作必须先制定方案，向用户说明后再执行
