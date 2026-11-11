# Transversity / Collins fitting workflow

采用 **Julia 数值内核 + Python/Snakemake 完整任务调度**（ADR-001）。目标是联合 SIDIS、SIA、STAR hadron-in-jet 与可选 lattice tensor charge，完成历史复现、数据更新、联合拟合、演化和可复现产物。

当前为 **可执行合成闭环，生产物理拟合尚未开放**。66 个合成测试点检验协方差/nuisance、预测、AD、多起点优化、replicas、holdouts、缓存和 receipts；它们不是实验测量。KPSY15 作者原始程序与修改的 HOPPET 尚未找到。

## 运行

使用 Linux/WSL2、Python 3.12、Julia 1.12.6。环境和运行目录放在 Linux 文件系统。默认一任务一进程，Julia/BLAS/OpenMP 线程数均为 1。

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
julia --project=numerics -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
snakemake --snakefile workflow/Snakefile --cores 2 --config mode=closure run_dir=runs/closure
```

输出包括 receipts、数值验证报告、fit/replica 数组、逐 bin CSV、合成闭环图和 `products/summary.json`。失败/BLOCKED 停止下游，`failures/` 保留证据。新版本使用新运行目录。也可用 `collinsflow run --mode closure --run-dir runs/local` 顺序执行同一 DAG。

```bash
python -m pytest -q
python scripts/plan_workflow.py validate
julia --project=numerics numerics/test/runtests.jl
make -C fortran
julia --project=numerics numerics/bin/check_bridge.jl "$PWD/fortran/libcollins_operator.so"
collinsflow readiness  # 当前预期退出 2，并列出生产缺项
```

CI 在 Linux 运行测试、Fortran ABI 和完整 Snakemake 闭环，并保留 30 天 artifacts。`main` 使用 PR + squash merge；见 [Git 规则](docs/git-policy.md) 和 [服务端实际配置](config/github-state.json)。

## 进度和科学边界

以 [Issues](https://github.com/Renoir-J/transversity-collins-fit/issues)、[Milestones](https://github.com/Renoir-J/transversity-collins-fit/milestones) 和 [实现状态](docs/implementation-status.md) 管理进度。生产入口受 [physics contract](config/physics-contract.json) 限制，不能把 synthetic receipt 或一个布尔开关当成科学验收。

- [Julia/Fortran 职责与验证](docs/numerical-workflow.md)
- [数值 artifact 与数据接口](docs/kernel-contract.md)
- [12 阶段 DAG 与 7 个研究场景](config/workflow.json)
- [数据候选、替换和重叠策略](config/data-registry.json)
- [公开代码审查](references/code-review.json)和[八个仓库的固定版本](references/source-inventory.json)
- [完整理论与文献研究 Page](https://chatgpt.com/space/page_b154ab4fd11c819183553e02ea97830b)

KPSY CSS NLL′、JAM Gaussian/log-log、固定 nf 的 LO moment 检查和 FeynFacet NLO hard functions 属于不同层级/约定，不能互相冒充复现。`external/`、原始下载、私信备份、环境及 runs 不进入 git。见 [第三方许可说明](THIRD_PARTY_NOTICES.md)。
