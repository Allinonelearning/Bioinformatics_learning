#!/bin/bash

# 定义源文件夹和目标文件夹
SOURCE_DIR="/home/us108/Knl1/all_pep"
DEST_DIR="/home/us108/Knl1/NASP"

# 检查目标文件夹是否存在，如果不存在则创建
if [ ! -d "$DEST_DIR" ]; then
  mkdir -p "$DEST_DIR"
fi

# 遍历源文件夹中的所有.pep文件
for FILE in "$SOURCE_DIR"/*.pep; do
    ../hmmer/hmmer-3.3.2/hmmsearch --domtblout "$FILE".out PF10516.hmm $FILE

     # 检查.out文件是否存在
    if [ -f "${FILE}.out" ]; then
        # 将生成的.out文件移动到目标文件夹
        mv "${FILE}.out" "$DEST_DIR"
    else
        echo "Failed to generate .out file for $FILE"
    fi

done

echo "All .pep files have been processed and moved to $DEST_DIR"
