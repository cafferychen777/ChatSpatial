#!/bin/bash

# LaTeX Workflow Automation Script
# 自动化 LaTeX 文档编译、预览和同步

set -e

# 配置
MAIN_TEX=""
BUILD_DIR="build"
WATCH_MODE=false

# 颜色
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 函数：查找主 TeX 文件
find_main_tex() {
    # 查找包含 \documentclass 的文件
    MAIN_FILES=$(grep -l "\\\\documentclass" *.tex 2>/dev/null || true)
    
    if [ -z "$MAIN_FILES" ]; then
        echo -e "${RED}错误: 未找到主 TeX 文件${NC}"
        exit 1
    fi
    
    # 如果只有一个文件，自动选择
    FILE_COUNT=$(echo "$MAIN_FILES" | wc -l)
    if [ "$FILE_COUNT" -eq 1 ]; then
        MAIN_TEX=$MAIN_FILES
    else
        echo -e "${YELLOW}找到多个 TeX 文件：${NC}"
        select tex_file in $MAIN_FILES; do
            MAIN_TEX=$tex_file
            break
        done
    fi
    
    echo -e "${GREEN}使用主文件: $MAIN_TEX${NC}"
}

# 函数：编译 LaTeX
compile_latex() {
    local tex_file=${1:-$MAIN_TEX}
    
    if [ -z "$tex_file" ]; then
        find_main_tex
        tex_file=$MAIN_TEX
    fi
    
    echo -e "${BLUE}编译 $tex_file ...${NC}"
    
    # 创建构建目录
    mkdir -p $BUILD_DIR
    
    # 第一次编译
    pdflatex -output-directory=$BUILD_DIR -interaction=nonstopmode "$tex_file" || {
        echo -e "${RED}PDFLaTeX 编译失败${NC}"
        return 1
    }
    
    # 检查是否有参考文献
    if grep -q "\\\\bibliography\\|\\\\addbibresource" "$tex_file"; then
        echo -e "${BLUE}处理参考文献...${NC}"
        cd $BUILD_DIR
        bibtex "${tex_file%.tex}" || echo -e "${YELLOW}BibTeX 警告（可能正常）${NC}"
        cd ..
        
        # 需要两次额外编译以更新引用
        pdflatex -output-directory=$BUILD_DIR -interaction=nonstopmode "$tex_file"
        pdflatex -output-directory=$BUILD_DIR -interaction=nonstopmode "$tex_file"
    fi
    
    echo -e "${GREEN}✓ 编译完成: $BUILD_DIR/${tex_file%.tex}.pdf${NC}"
}

# 函数：清理临时文件
clean_latex() {
    echo -e "${YELLOW}清理临时文件...${NC}"
    rm -f *.aux *.log *.out *.toc *.bbl *.blg *.synctex.gz
    rm -rf $BUILD_DIR
    echo -e "${GREEN}✓ 清理完成${NC}"
}

# 函数：监视模式
watch_latex() {
    if ! command -v fswatch &> /dev/null; then
        echo -e "${RED}错误: fswatch 未安装${NC}"
        echo "请运行: brew install fswatch"
        exit 1
    fi
    
    find_main_tex
    echo -e "${BLUE}监视模式启动 - 监控 *.tex 文件变化${NC}"
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}"
    
    # 初始编译
    compile_latex
    
    # 监视文件变化
    fswatch -o *.tex | while read f; do
        echo -e "${YELLOW}检测到文件变化，重新编译...${NC}"
        compile_latex
    done
}

# 函数：打开 PDF
open_pdf() {
    find_main_tex
    PDF_FILE="$BUILD_DIR/${MAIN_TEX%.tex}.pdf"
    
    if [ ! -f "$PDF_FILE" ]; then
        echo -e "${YELLOW}PDF 不存在，先编译...${NC}"
        compile_latex
    fi
    
    echo -e "${BLUE}打开 PDF: $PDF_FILE${NC}"
    
    if [[ "$OSTYPE" == "darwin"* ]]; then
        open "$PDF_FILE"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        xdg-open "$PDF_FILE"
    else
        echo -e "${RED}不支持的操作系统${NC}"
    fi
}

# 函数：Git 同步
sync_git() {
    echo -e "${BLUE}同步 Git 仓库...${NC}"
    
    # 检查是否在 Git 仓库中
    if ! git rev-parse --git-dir > /dev/null 2>&1; then
        echo -e "${RED}错误: 不在 Git 仓库中${NC}"
        return 1
    fi
    
    # 拉取最新更改
    echo -e "${YELLOW}拉取远程更改...${NC}"
    git pull origin main || git pull origin master
    
    # 添加所有更改
    git add -A
    
    # 检查是否有更改
    if git diff --staged --quiet; then
        echo -e "${GREEN}✓ 没有本地更改需要提交${NC}"
    else
        # 提交更改
        read -p "输入提交信息: " commit_msg
        git commit -m "$commit_msg"
        
        # 推送到远程
        echo -e "${YELLOW}推送到远程仓库...${NC}"
        git push origin main || git push origin master
        echo -e "${GREEN}✓ 同步完成${NC}"
    fi
}

# 函数：创建新文档
new_document() {
    read -p "输入文档名称（不含.tex）: " doc_name
    
    cat > "${doc_name}.tex" << 'EOF'
\documentclass[12pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage[margin=1in]{geometry}

\title{Document Title}
\author{Your Name}
\date{\today}

\begin{document}

\maketitle

\begin{abstract}
This is the abstract of your document.
\end{abstract}

\section{Introduction}
Your content starts here.

\section{Methodology}

\section{Results}

\section{Conclusion}

\bibliographystyle{plain}
% \bibliography{references}

\end{document}
EOF
    
    echo -e "${GREEN}✓ 创建文档: ${doc_name}.tex${NC}"
}

# 函数：显示帮助
show_help() {
    cat << EOF
LaTeX 工作流自动化脚本

用法: $0 [选项]

选项:
    compile     编译 LaTeX 文档
    clean       清理临时文件
    watch       监视模式（自动重新编译）
    open        打开生成的 PDF
    sync        同步 Git 仓库
    new         创建新文档模板
    help        显示此帮助信息

示例:
    $0 compile  # 编译当前目录的 LaTeX 文档
    $0 watch    # 启动监视模式
    $0 sync     # 同步到 Overleaf/GitHub
EOF
}

# 主程序
case "${1:-help}" in
    compile) compile_latex ;;
    clean) clean_latex ;;
    watch) watch_latex ;;
    open) open_pdf ;;
    sync) sync_git ;;
    new) new_document ;;
    help) show_help ;;
    *) show_help; exit 1 ;;
esac