# QBO 实验

[English](README.md) · [架构说明](docs/architecture.zh-CN.md) · [数据与运行流程](docs/workflows.zh-CN.md)

这个项目使用 CESM/WACCM，研究大气对准两年振荡（QBO）的响应如何随气候背景变化。实验将同一组西风、东风 QBO 目标分别施加到早期和晚期海温、海冰背景中，比较两种背景下的响应。

仓库包含输入场构建、对流层顶边界权重、强迫诊断，以及 CAM 接入所用的 Fortran 验证程序。

## 实验设计

| | 早期背景：1981–1994 | 晚期背景：1995–2024 |
|---|---|---|
| 西风 QBO | W × early | W × late |
| 东风 QBO | E × early | E × late |

QBO 目标来自 MERRA-2 赤道月平均纬向风。按 1981–2024 年各年 3 月 70 hPa 风速距平排序，取最高和最低各 8 年，构建 W、E 两组月合成。每组包含 7 个气压层、13 个月节点，覆盖前一年 9 月至当年 9 月。同一相位在两种背景中使用相同的目标文件。

海温和海冰背景来自 HadISST 月气候态，经周期经度最近邻插值映射到基准 FV 网格；源数据缺测处采用对应位置的基准场。

## 安装与测试

需要 Python 3.10 或以上版本：

```bash
python -m venv .venv
```

Linux/macOS 使用 `source .venv/bin/activate` 激活环境，PowerShell 使用 `.venv\Scripts\Activate.ps1`，然后在仓库根目录执行：

```bash
python -m pip install -e ".[test]"
python -m pytest
```

Python 测试在本地生成小型数组和 NetCDF 文件，覆盖输入维度、经度周期拼接、边界权重、强迫诊断、纬向与非纬向分解，以及 CAM 源码准备。

## 构建输入

将 MERRA-2 月风速表、HadISST 压缩数据和 CESM 海温基准文件放到 `data/`。所需变量和处理步骤见[数据与运行流程](docs/workflows.zh-CN.md)。

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

构建程序输出 NetCDF 文件及带输入、输出哈希的 JSON 元数据。QBO 文件采用 `qbo(time, level)` 存储顺序，对应 WACCM Fortran 端的 `u_inp(level, time)` 数组。

## 使用边界算子

```python
import numpy as np
from qbo.boundary import safe_pressure, column_taper
from qbo.zonal import zonal_envelope

# 三组对流层顶诊断、四个经度点，气压单位为 hPa。
pressure = np.array([[95., 100., 105., 110.]] * 3)
safe = safe_pressure(pressure, np.ones_like(pressure))
column_weights = column_taper(80., safe)
ring_weight = zonal_envelope(column_weights)
```

逐柱算子以三组有效对流层顶中最严格的边界为准，计算平滑衰减权重。纬圈算子再取所有经度的最小值，使同一纬圈的强迫权重一致，并满足每一柱的边界约束。

## 代码结构

```text
src/qbo/          输入构建、数值算子与诊断
scripts/          CAM 源码准备和历史运行分析
native/           Fortran 边界、纬圈、NetCDF 与 MPI 契约测试
tests/            Python 回归测试
docs/             架构、运行流程与代码来源
vendor/cam/       固定版本的 CAM QBO 源码及原许可
```

Fortran 检查需要 GNU Fortran 和 MPI 工具链：

```bash
python scripts/check_native.py --output outputs/native-check
```

该命令验证 Python 与 Fortran 边界计算的一致性、不同 MPI 分配下的纬圈归约，以及 12 个时间步和逐柱覆盖测试案例。

## 当前进度

输入构建、离线算子、带诊断的 v4 CAM 适配器，以及独立的 v5 MPI 契约已有测试覆盖。v5 的 prepare/reduce/resume 接入方案见[架构说明](docs/architecture.zh-CN.md#cam-接入)。完整模式中的背景响应实验是下一阶段工作。

历史分析脚本通过 `QBO_WORKSPACE` 读取实验档案。大型输入数据、模式历史场、作业日志和本机连接配置保存在 Git 仓库之外。

固定版本的 CAM 源码保留其[原许可](vendor/cam/LICENSE.txt)。
