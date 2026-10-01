# QBO 气候背景对照实验

[English](README.md) · [架构说明](docs/architecture.zh-CN.md) · [数据与运行流程](docs/workflows.zh-CN.md)

用于 CESM/WACCM 准两年振荡（QBO）实验的输入构建与强迫诊断代码。实验固定西风和东风两组 QBO 目标，分别搭配早期、晚期海温和海冰背景，研究大气响应随气候背景的变化。

## 实验设计

| | 早期背景：1981–1994 | 晚期背景：1995–2024 |
|---|---|---|
| 西风 QBO | W × early | W × late |
| 东风 QBO | E × early | E × late |

QBO 目标由 MERRA-2 赤道月平均纬向风构建。按 1981–2024 年各年 3 月 70 hPa 风速距平排序，取最高和最低各 8 年，分别合成 W、E 两组目标。每组包含 7 个气压层、13 个月节点，从入选年的前一年 9 月延续至当年 9 月。两种背景使用相同的相位目标文件。

海温和海冰背景采用两个时期的 HadISST 月气候态，通过最近邻插值映射到基准 FV 网格，经度按周期处理。源数据缺测处使用对应位置的基准场填补。

## 安装

需要 Python 3.10 或以上版本：

```bash
git clone https://github.com/ARETE-zzwl/qbo.git
cd qbo
python -m venv .venv
```

Linux/macOS 使用 `source .venv/bin/activate` 激活环境，PowerShell 使用 `.venv\Scripts\Activate.ps1`，然后在仓库根目录执行：

```bash
python -m pip install -e ".[test]"
```

## 运行示例

```bash
python -m pip install -e ".[test,demo]"
qbo-demo --output outputs/demo
```

示例在本地生成风场、海温、海冰和对流层顶合成数据，运行两类输入构建程序，输出目标文件检查、边界权重，以及 PNG、SVG、PDF 格式的总览图。打开 `outputs/demo/README.md` 查看结果。每次运行使用新的输出目录。

分析自己的快照文件：

```bash
qbo-diagnose --input data/snapshot.nc --output outputs/diagnostics
```

[快照格式](docs/workflows.zh-CN.md#诊断单个快照)要求提供完整经度网格上的气压场。[科研方向](docs/research-directions.zh-CN.md)整理了下一步模式验证和可深入的问题。

## 构建输入

将 MERRA-2 月风速表、HadISST 压缩数据和 CESM 海温基准文件放到 `data/`。所需变量和处理步骤见[数据与运行流程](docs/workflows.zh-CN.md)。

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

输出包括 NetCDF 文件和 JSON 元数据。QBO 文件采用 `qbo(time, level)` 存储顺序，与 WACCM 的 Fortran 读取接口对应。

## 边界权重

逐柱权重取三组对流层顶诊断中最严格的边界，在层底接近该边界时平滑衰减至零。沿经度取最小值后，得到同一纬圈共用的权重。

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

## 测试

```bash
python -m pytest
```

测试覆盖输入构建、NetCDF 维度、经度映射、边界权重、强迫诊断和 CAM 源码准备，所需测试数据在本地生成。

[GitHub Actions](https://github.com/ARETE-zzwl/qbo/actions/workflows/tests.yml) 在 Linux、Windows 和 Python 3.10、3.13 上运行测试与示例，并在 Linux 上编译、运行 Fortran/MPI 检查。每次运行保留报告和编译日志。

安装 GNU Fortran 和 MPI 后，可运行原生测试：

```bash
python scripts/check_native.py --output outputs/native-check
```

程序对照 Fortran 内核与 Python 参考值，并检查纬圈归约和 MPI 逐柱覆盖。每次运行使用新的输出目录。

## 代码结构

```text
src/qbo/          输入构建、数值算子与诊断
scripts/          CAM 源码准备和历史运行分析
native/           Fortran 边界、纬圈、NetCDF 与 MPI 契约测试
tests/            Python 回归测试
docs/             架构、运行流程与代码来源
vendor/cam/       固定版本的 CAM QBO 源码及原许可
```

v5 的 CAM 驱动接入目前处于设计阶段，方案见 [CAM 接入](docs/architecture.zh-CN.md#cam-接入)。分析已有模式运行时，通过 `QBO_WORKSPACE` 指定实验档案目录，具体用法见[数据与运行流程](docs/workflows.zh-CN.md#分析已有实验档案)。

`vendor/cam/` 保存适配器使用的 CAM QBO 源码及其[原许可](vendor/cam/LICENSE.txt)。
