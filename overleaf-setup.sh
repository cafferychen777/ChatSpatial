#!/bin/bash

# Overleaf Git Integration Setup Script
# 用于配置和管理 Overleaf 项目的本地 Git 同步

set -e

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}=== Overleaf Git Integration Setup ===${NC}"

# 函数：克隆 Overleaf 项目
clone_overleaf() {
    echo -e "${YELLOW}克隆 Overleaf 项目到本地${NC}"
    read -p "请输入 Overleaf 项目的 Git URL (在项目菜单中找到): " OVERLEAF_URL
    read -p "请输入本地项目文件夹名称: " PROJECT_NAME
    
    if [ -d "$PROJECT_NAME" ]; then
        echo -e "${RED}错误: 文件夹 $PROJECT_NAME 已存在${NC}"
        exit 1
    fi
    
    git clone "$OVERLEAF_URL" "$PROJECT_NAME"
    echo -e "${GREEN}✓ 项目克隆成功${NC}"
    
    cd "$PROJECT_NAME"
    
    # 设置 Git 配置
    git config core.autocrlf input
    git config pull.rebase false
    
    echo -e "${GREEN}✓ Git 配置完成${NC}"
}

# 函数：配置 GitHub 远程仓库（双远程）
setup_github_remote() {
    echo -e "${YELLOW}设置 GitHub 作为第二个远程仓库${NC}"
    read -p "是否要添加 GitHub 远程仓库？(y/n): " ADD_GITHUB
    
    if [ "$ADD_GITHUB" = "y" ]; then
        read -p "请输入 GitHub 仓库 URL: " GITHUB_URL
        git remote add github "$GITHUB_URL"
        echo -e "${GREEN}✓ GitHub 远程仓库添加成功${NC}"
        
        # 显示所有远程仓库
        echo -e "${YELLOW}当前远程仓库：${NC}"
        git remote -v
    fi
}

# 函数：同步到 Overleaf
push_to_overleaf() {
    echo -e "${YELLOW}推送更改到 Overleaf${NC}"
    git add -A
    read -p "请输入提交信息: " COMMIT_MSG
    git commit -m "$COMMIT_MSG"
    git push origin main
    echo -e "${GREEN}✓ 已推送到 Overleaf${NC}"
}

# 函数：从 Overleaf 拉取
pull_from_overleaf() {
    echo -e "${YELLOW}从 Overleaf 拉取最新更改${NC}"
    git pull origin main
    echo -e "${GREEN}✓ 已从 Overleaf 拉取最新版本${NC}"
}

# 函数：创建 LaTeX 编译环境检查
check_latex_env() {
    echo -e "${YELLOW}检查 LaTeX 环境${NC}"
    
    if command -v pdflatex &> /dev/null; then
        echo -e "${GREEN}✓ pdflatex 已安装${NC}"
    else
        echo -e "${RED}✗ pdflatex 未安装${NC}"
        echo "  建议安装 MacTeX (Mac) 或 TeX Live (Linux)"
    fi
    
    if command -v bibtex &> /dev/null; then
        echo -e "${GREEN}✓ bibtex 已安装${NC}"
    else
        echo -e "${RED}✗ bibtex 未安装${NC}"
    fi
}

# 主菜单
main_menu() {
    echo ""
    echo "请选择操作："
    echo "1) 克隆新的 Overleaf 项目"
    echo "2) 设置 GitHub 远程仓库"
    echo "3) 推送更改到 Overleaf"
    echo "4) 从 Overleaf 拉取更改"
    echo "5) 检查 LaTeX 环境"
    echo "6) 退出"
    
    read -p "选择 (1-6): " CHOICE
    
    case $CHOICE in
        1) clone_overleaf ;;
        2) setup_github_remote ;;
        3) push_to_overleaf ;;
        4) pull_from_overleaf ;;
        5) check_latex_env ;;
        6) echo -e "${GREEN}再见！${NC}"; exit 0 ;;
        *) echo -e "${RED}无效选择${NC}"; main_menu ;;
    esac
    
    main_menu
}

# 运行主菜单
main_menu