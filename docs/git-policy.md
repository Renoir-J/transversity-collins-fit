# Git 与研究版本规则

仓库默认私有。`main` 保存可复现、通过 CI 的版本；用 `feat/`、`fix/`、`theory/`、`data/`、`docs/` 分支提交 PR，squash merge 后删除分支。每项研究工作以 Issue 描述验收证据并归入 milestone。单人仓库不强制一个本人无法提供的额外 review；涉及理论约定的变更仍需明确记录科学审查状态。

本地安装规则：`git config --local core.hooksPath .githooks`，`git config --local pull.ff only`，`git config --local fetch.prune true`。pre-commit 检查 staged 内容；pre-push 拒绝直接推 main 并运行无需额外依赖的规格/内容检查。完整数值测试由 Linux CI 执行。hook 可被本地操作者跳过，不能等同于服务端保护。

目标服务端规则：PR 合并、CI `validation` 通过、讨论解决、线性历史、禁止 force push / 删除 main，并适用于管理员。GitHub Free 的私有仓库可能无法启用这类规则；实际启用结果记录在 `config/github-state.json`。不会为获得规则而把研究仓库改成公开，也不会自动购买升级。[GitHub protected branches](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches)

提交代码、锁文件、小型合成测试与来源索引。下载的第三方代码、原始数据、运行结果、环境、Slack 导出和凭据不入 git。保留上游许可证；新代码尚未选择公开许可证。运行产物保存到 CI artifacts 或独立归档，附代码、环境、输入与输出 hash。每次拟合使用新目录，receipt 不能从其他 scope 或旧版本复制。

数据替换必须删除旧 likelihood。新增数据、符号、尺度、jet R/axis、积分、模型及演化变更须触发下游失效。发布物理结果前，所有 production gates 必须通过；合成闭环通过只说明所测软件流程可运行。

