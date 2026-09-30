# Source organization / 代码整理对应表

The numerical functions below were moved without changing their calculations. Input builders now share one SHA-256 helper, and tests import the installed package.

下列数值函数保留原有计算方式，改按职责组织。输入构建程序共用 SHA-256 函数，测试直接导入安装后的包。

| Original file / 原文件 | Current location / 当前位置 |
|---|---|
| `build_same_qbo_targets.py` | `src/qbo/targets.py` |
| `build_early_late_hadisst_backgrounds.py` | `src/qbo/backgrounds.py` |
| `validate_same_qbo_targets.py` | `src/qbo/validate.py` |
| `scripts/qbo_boundary_v4.py` | `src/qbo/boundary.py` |
| `scripts/qbo_zonal_v5.py` | `src/qbo/zonal.py` |
| Numerical helpers in `audit_remote_qbo.py` | `src/qbo/forcing.py` |
| `retention_metrics`, `operator_metrics` | `src/qbo/metrics.py` |
| Namelist and cutoff helpers in `audit_v3_preconditions.py` | `src/qbo/preconditions.py` |
| `scripts/native_boundary_v4.py` | `src/qbo/native.py` |
| `scripts/prepare_cam_adapter_v4.py` | `scripts/prepare_cam_adapter.py` |
| Archive analysis entry points | `scripts/analysis/` |
| Boundary and CAM adapter Fortran sources | `native/boundary/`, `native/cam_adapter/` |
| Zonal R2 driver and isolated coverage contract | `native/zonal/`, `native/contract/` |

The pinned CAM baseline and native numerical sources retain their original bytes. The source map records hashes before Python import/path edits. Site-specific SSH collectors, Slurm launchers, monitoring state and experiment logs remain in the original workspace.

固定版本的 CAM 基准和原生数值源码保留原始字节。来源表中的哈希记录 Python 导入和路径调整前的版本。面向原服务器的 SSH 收集脚本、Slurm 提交脚本、监控状态和实验日志保留在原工作目录。
