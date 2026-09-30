# 数据与运行流程

[English](workflows.md) · [README](../README.zh-CN.md)

安装 `.[test]` 后，在仓库根目录执行以下命令。

## 输入文件

| 文件 | 所需内容 |
|---|---|
| `QBO_MERRA2-Uvals_00N_GSFC.txt` | ASCII 月数据表，包含 `P (hPa):` 气压头和 `YYYYMM` 数据行；选用 10、20、30、40、50、70、100 hPa |
| `HadISST_sst.nc.gz` | `sst(time, latitude, longitude)` 和 CF 时间坐标 |
| `HadISST_ice.nc.gz` | 相同网格上的 `sic(time, latitude, longitude)`，海冰覆盖率范围为 0–1 |
| `sst_baseline.nc` | `lat`、`lon`、`time`、`date`、`datesec`、`SST_cpl`、`ice_cov`；两个场均需 `long_name` 和 `units` 属性 |

QBO 年份筛选范围为 1981–2024。合成从入选年 3 月之前的 9 月开始，因此输入还需覆盖前一年的 9–12 月。HadISST 需覆盖 1981 年 1 月至 2024 年 12 月的全部月份。海温基准文件提供目标网格上的 12 个月记录。

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

目标构建程序输出 `qbo_target_W.nc`、`qbo_target_E.nc`、`selection.csv` 和 `target_metadata.json`。NetCDF 日期采用模式第 1、2 年的月初节点，末端再次包含 9 月。因此，3 月节点是时间插值节点，其值与整月积分平均值的含义不同。

背景构建程序输出早期、晚期两份气候态及 `background_metadata.json`。海冰覆盖率截断到 0–1，主字段和 `_prediddle` 字段初始值相同。CESM 的月平均保持调整属于后续模式输入处理步骤。

## 原生测试

Linux 环境的 `PATH` 中有 `gfortran`、`mpifort`、`mpiexec` 后，执行：

```bash
python scripts/check_native.py --output outputs/native-check
```

每次运行使用新的输出目录。程序保存编译与运行日志、测试向量、数值对照和 JSON 汇总，依次检查：

1. Fortran 逐柱内核与 10,047 组 Python 向量的一致性，包含边界值和无效输入。
2. 小型纬圈数据在单/双 rank、cyclic、block、single-owner 分配下的结果，以及覆盖不完整的案例。
3. `native/contract/cases.json` 定义的 12 个双 rank 契约案例。

`native/netcdf/read_qbo_native.f90` 用于具备 NetCDF Fortran 库的环境。Python 测试也会调用 NetCDF C API 验证目标文件的维度顺序。

## 准备 CAM 源码

```bash
python scripts/prepare_cam_adapter.py outputs/cam-adapter
```

命令创建新目录，包含 `SourceMods/src.cam/`、原生测试源码、`manifest.json` 和 `SHA256SUMS`。基准为 `vendor/cam/qbo.F90`，对应归档的 CAM `cam_cesm2_1_rel_60` 源码。文件准备在本地完成，CESM case 的编译与运行使用所在计算环境的工具链和配置。

## 分析已有实验档案

分析脚本保留了 2026 年 9 月的档案布局和输入哈希。将 `QBO_WORKSPACE` 指向包含原 `outputs/`、`scripts/`、实验协议和 `source_paths.json` 的目录：

```bash
export QBO_WORKSPACE=/path/to/qbo-archive
python scripts/analysis/evaluate_boundary_v4.py --output outputs/column-assessment.json
python scripts/analysis/audit_qbo_operator_scope.py --output outputs/operator-assessment.json
python scripts/analysis/assess_zonal_v5.py --output outputs/zonal-assessment.json
```

PowerShell 中使用 `$env:QBO_WORKSPACE = 'D:/path/to/qbo-archive'`。档案内的 `source_paths.json` 指向原始输入实验和 HadISST 目录，安装、输入构建及回归测试均无需该文件。

| 脚本 | 读取的档案内容 |
|---|---|
| `audit_inputs.py` | 原始/修正目标和背景场、固定版本的 CAM 源码包 |
| `audit_remote_qbo.py` | 原生读取输出、垂直网格和作业检查点 |
| `audit_qbo_technical.py` | 短积分历史场和耦合器结束日志 |
| `audit_v3_preconditions.py` | case 配置、namelist 和编译记录 |
| `screen_qbo_boundaries.py` | 用于硬边界比较的历史场 |
| `evaluate_boundary_v4.py` | 归档状态及 v3 分析结果 |
| `audit_qbo_operator_scope.py` | v4 分析结果和源码指纹 |
| `assess_zonal_v5.py` | v4 分析结果及同一组 6 份归档状态 |
| `summarize_remote_increment.py` | 运行报告和 Slurm 记账 |

带 `--output` 参数的脚本可指定新输出位置。输入、原生读取、技术积分及汇总审计会将报告写回既定档案位置；需要保留原报告时，使用档案副本。
