# GitHub Actions 优化指南

## 🚨 当前问题
- GitHub Actions 免费额度已用完（3000/3000 分钟）
- ChatSpatial 消耗了 $22.64/月
- 需要等待下月重置或升级付费计划

## ✅ 已实施的优化

### 1. 整合 Workflow（减少 85% 运行）
**之前：** 7 个 workflow，每次 push 都运行
**现在：** 1 个主 workflow + 其他手动触发

### 2. 优化后的结构
```
main.yml          - 主 CI（自动运行，轻量级）
test.yml          - 完整测试（手动触发）
tests.yml         - 单元测试（手动触发）
ci.yml            - 详细 CI（手动触发）
minimal-test.yml  - 最小测试（手动触发）
quality.yml       - 代码质量（手动触发）
release.yml       - 发布（标签触发）
```

### 3. 节省措施
- ❌ 移除每日定时任务（节省 30×30 = 900 分钟/月）
- ❌ 移除多 Python 版本测试（节省 50%）
- ✅ 路径过滤（忽略文档更改）
- ✅ 缓存依赖（减少安装时间）

## 📊 预期效果
- **之前：** ~3000+ 分钟/月
- **之后：** ~200-300 分钟/月
- **节省：** ~90% Actions 使用量

## 🔧 手动运行测试
```bash
# 在 GitHub Actions 页面手动触发
gh workflow run test.yml
gh workflow run tests.yml
gh workflow run quality.yml
```

## 💡 其他建议

### 短期解决方案
1. 等待 9月1日 额度重置
2. 或升级到 GitHub Team（3000 → 3000 分钟）

### 长期解决方案
1. 使用自托管 runner（免费无限）
2. 只在 PR 合并前运行完整测试
3. 使用 act 本地运行 Actions

## 📈 监控使用量
访问 [GitHub Billing](https://github.com/settings/billing) 查看使用情况

---
更新时间：2025-08-26