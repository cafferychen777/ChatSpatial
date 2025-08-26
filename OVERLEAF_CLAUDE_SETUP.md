# Claude Code + Overleaf 集成指南

## 🚀 快速开始

虽然 Claude Code 没有直接的 Overleaf 插件，但通过 Git 可以实现完美的工作流集成。

### 前置条件
- ✅ Git 已配置（用户名：Chen Yang，邮箱：cafferychen777@tamu.edu）
- ✅ SSH 密钥已设置（id_ed25519, id_rsa）
- ✅ Claude Code 已安装

## 📦 安装的工具

我已经为你创建了三个强大的脚本工具：

### 1. `overleaf-setup.sh` - Overleaf 项目管理
```bash
./overleaf-setup.sh
```
功能：
- 克隆 Overleaf 项目到本地
- 设置 GitHub 作为第二远程仓库
- 推送/拉取 Overleaf 更改
- 检查 LaTeX 环境

### 2. `latex-workflow.sh` - LaTeX 自动化工作流
```bash
./latex-workflow.sh [命令]
```
命令：
- `compile` - 编译 LaTeX 文档
- `watch` - 监视模式（文件改动自动重编译）
- `clean` - 清理临时文件
- `open` - 打开生成的 PDF
- `sync` - 同步到 Git 仓库
- `new` - 创建新文档模板

### 3. `claude-latex-helper.sh` - Claude Code LaTeX 助手
```bash
./claude-latex-helper.sh
```
功能：
- LaTeX 公式助手
- 表格生成器
- 错误诊断工具
- 常用包参考
- 快速 Overleaf 同步
- 生成 Claude 提示词

## 🔄 推荐工作流程

### 方案一：Overleaf Git 集成（推荐）

1. **获取 Overleaf 项目的 Git URL**
   - 打开 Overleaf 项目
   - 点击菜单 → Git
   - 复制 Git URL

2. **克隆到本地**
   ```bash
   ./overleaf-setup.sh
   # 选择 1 克隆项目
   ```

3. **在 Claude Code 中编辑**
   ```bash
   # 让 Claude 帮你写 LaTeX
   claude "帮我写一个关于机器学习的介绍章节"
   
   # 或者修复错误
   claude "修复这个 LaTeX 编译错误"
   ```

4. **本地编译测试**
   ```bash
   ./latex-workflow.sh compile
   ```

5. **同步回 Overleaf**
   ```bash
   ./latex-workflow.sh sync
   ```

### 方案二：监视模式开发

1. **启动监视模式**
   ```bash
   ./latex-workflow.sh watch
   ```

2. **在另一个终端使用 Claude Code**
   ```bash
   claude "添加一个算法伪代码到 methods.tex"
   ```

3. **文件保存后自动编译并预览**

### 方案三：使用 LaTeX 助手

1. **运行助手**
   ```bash
   ./claude-latex-helper.sh
   ```

2. **选择需要的功能**
   - 生成复杂公式
   - 创建表格
   - 诊断错误

## 💡 Claude Code 使用技巧

### LaTeX 相关提示词示例

```bash
# 写公式
claude "写一个贝叶斯公式的 LaTeX 代码"

# 创建表格
claude "创建一个 3x4 的实验结果对比表格"

# 修复错误
claude "修复 Undefined control sequence 错误"

# 优化文档
claude "优化这个 LaTeX 文档的结构和格式"

# 添加引用
claude "帮我设置 BibTeX 并添加引用"
```

### 高级工作流

1. **批量处理**
   ```bash
   claude "将所有的 equation 环境改为 align 环境"
   ```

2. **自动格式化**
   ```bash
   claude "格式化这个 LaTeX 文件，确保一致的缩进和换行"
   ```

3. **生成图表**
   ```bash
   claude "用 TikZ 画一个神经网络架构图"
   ```

## ⚠️ 注意事项

1. **Overleaf 限制**
   - Git 集成是付费功能
   - 不支持分支操作
   - 同步需要手动触发

2. **本地环境**
   - 需要安装 LaTeX 发行版（MacTeX 或 TeX Live）
   - 推荐安装 `fswatch` 用于文件监视：`brew install fswatch`

3. **最佳实践**
   - 定期提交和同步
   - 使用有意义的提交信息
   - 在本地测试后再推送到 Overleaf

## 🛠️ 故障排除

### 问题：Overleaf 同步失败
```bash
# 检查远程仓库
git remote -v

# 重新设置 Overleaf 远程
git remote set-url origin [Overleaf Git URL]
```

### 问题：LaTeX 编译错误
```bash
# 使用诊断工具
./claude-latex-helper.sh
# 选择 3 - 错误诊断
```

### 问题：PDF 无法预览
```bash
# 清理并重新编译
./latex-workflow.sh clean
./latex-workflow.sh compile
./latex-workflow.sh open
```

## 📚 额外资源

- [Overleaf Git 集成文档](https://www.overleaf.com/learn/how-to/Git_integration)
- [Claude Code 官方文档](https://claude.ai/docs)
- [LaTeX 参考手册](https://en.wikibooks.org/wiki/LaTeX)

---

现在你可以开始使用了！运行 `./overleaf-setup.sh` 开始设置你的第一个项目。