# QBO 气候背景对照实验

[English](README.md) · [架构](docs/architecture.zh-CN.md) · [数据格式](docs/workflows.zh-CN.md)

CESM/WACCM 准两年振荡（QBO）实验的输入构建与强迫诊断工具。

代码从 MERRA-2 赤道风场构建西风、东风两组 QBO 目标，从 HadISST 构建早期（1981–1994）和晚期（1995–2024）的海温、海冰背景。两种背景共用同一组 QBO 目标。Python 与 Fortran 程序负责计算逐柱边界权重及其在各纬圈上的最小值。

## 安装

需要 Python 3.10 或以上版本。

```bash
git clone https://github.com/ARETE-zzwl/qbo.git
cd qbo
python -m venv .venv
```

Linux/macOS 使用 `source .venv/bin/activate` 激活环境，PowerShell 使用 `.venv\Scripts\Activate.ps1`，然后安装：

```bash
python -m pip install -e ".[test,demo]"
```

## 运行示例

```bash
qbo-demo --output outputs/demo
```

示例生成合成输入，运行输入构建和边界诊断。结果汇总在 `outputs/demo/README.md`，包括 NetCDF 数据、CSV/JSON 统计和 PNG、SVG、PDF 格式的图。每次运行使用新的输出目录。

## 构建实验输入

将源文件放到 `data/`，执行：

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

QBO 目标按 1981–2024 年各年 3 月 70 hPa 风速距平排序，取最高和最低各 8 年合成 W/E 目标。每组包含 7 个气压层和 13 个月节点，从入选年的前一年 9 月延续至当年 9 月。文件按 `qbo(time, level)` 存储，供 WACCM 的 Fortran 接口读取。

HadISST 月气候态通过最近邻插值映射到基准 FV 网格，经度按周期处理，缺测网格使用基准值填补。所需变量和输出文件见[数据格式](docs/workflows.zh-CN.md)。

## 诊断快照

```bash
qbo-diagnose --input data/snapshot.nc --output outputs/diagnostics
```

柱权重和纬圈权重写入 `weights.nc`，统计结果保存为 CSV 和 JSON。气压变量与网格要求见[快照格式](docs/workflows.zh-CN.md#诊断单个快照)。

## 测试

```bash
python -m pytest
```

安装 GNU Fortran 和 MPI 后：

```bash
python scripts/check_native.py --output outputs/native-check
```

[CI](https://github.com/ARETE-zzwl/qbo/actions/workflows/tests.yml) 在 Linux、Windows 上运行 Python 测试和示例，并在 Linux 上运行 Fortran/MPI 检查。

## 目录

```text
src/qbo/       输入构建、算子和诊断
scripts/       CAM 源码准备与档案分析
native/        Fortran、MPI 内核与测试
tests/         Python 测试
docs/          架构与用法
vendor/cam/    CAM 源码及原许可
```

CAM 适配器的使用和驱动接入见[架构文档](docs/architecture.zh-CN.md#cam-接入)。仓库中的 CAM 源码保留[原许可](vendor/cam/LICENSE.txt)。
