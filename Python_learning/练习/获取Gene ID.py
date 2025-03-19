# 从原始文件中读取内容
with open('gene_result.txt', 'r') as file:
    data = file.readlines()

# 用于存储保留的行
desired_lines = []

# 遍历每一行数据
for line in data:
    # 去除行首和行尾的空白字符
    cleaned_line = line.strip()
    # 检查行是否以数字和点号开头，如果是，表示这是序号行
    if any(cleaned_line.startswith(str(i) + '.') for i in range(1, 1915)):
        # 去除序号和点号，只保留内容
        cleaned_line = cleaned_line.split('.', 1)[-1].strip()
        desired_lines.append(cleaned_line)

# 将desired_lines中的内容写入新的文本文件
with open('Aft.txt', 'w') as output_file:
    for line in desired_lines:
        output_file.write(line + '\n')
