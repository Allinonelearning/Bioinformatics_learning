# 多态性是指同一种行为具有多种表现形式
# 静态方法，可以使用对象访问也可以使用类访问
# class animal(object):
#     @staticmethod
#     def eat(name):
#         print(f"{name}") #不用self
#         print("狗会吃饭")
# animal.eat("xing")
# pe = animal()
# pe.eat("xing")

# 类方法 @classmethod 第一个参数必须是类对象，一般是els
# class people(object): 
#     name = "xingjinyin"
#     @classmethod
#     def eat(cls):
#         print(f"{cls.name}会吃饭")
# people.eat()
# xing = people()
# xing.eat()
    

# 文件操作
# 1.open():创建一个file对象，默认是只读模式
# 2.read(n)
# 3.write()
# 4.close()
# 文件属性
# 文件名.name
# 文件名.mode
# 文件名.closed


# 以下教程来自https://omicstutorials.com/python-strings-for-bioinformatics-from-basics-to-applications/
# 两个字符进行连接
# str1 = "hello"
# str2 = "xinngjinyin"
# Greeting = str1 + "" + str2
# print(Greeting)


# 获取字符串的长度
# string = "Python"
# length = len(string)
# print(length)

#索引 访问字符串中的单个字符
# string = "Python"
# print(string[0])

#切片
# substring = string[2:5]
# print(string[::2])
# print(substring)


#格式化 使用占位符format()方法格式化字符串
# name = "Alice"
# age = 30
# message = "My name is {} and I am {} years old.".format(name,age)
# print(message)

# 转换
# string = "hello World"
# print(string.lower())
# print(string.upper())

# 拆分,将字符串拆分成列表
# string = "apple,banana,peach"
# fruits = string.split(",")
# print(fruits)


# 连接
# fruits = ['apple', 'banana', 'peach']
# string = ",".join(fruits)
# print(type(string))

# Biopython 学习
# from Bio import SeqIO
# fast_file = "sequence.fasta"
# sequences = []
# sequence_id = "Sequence_1"
# sequence = "ATCGCGAAACGTACGATCGTACT"
# with open("learning.fasta","w") as file:
#     file.write(f">{sequence_id}\n")
#     file.write(f"{sequence}\n")

# sequences = {
# "Sequence_1": "ATCGATCGATCG",
# "Sequence_2": "GCTAGCTAGCTA",
# "Sequence_3": "TATGTATGTATG"
# }
# # Open the output file for writing
# with open("output.fasta", "w") as file:
# # Write each sequence in FASTA format
#     for sequence_id, sequence in sequences.items():
#         file.write(f">{sequence_id}\n") # Write the sequence ID
#         file.write(f"{sequence}\n") # Write the sequence






# 成对序列比对
# def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
# # Initialize the scoring matrix
#     rows = len(seq1) + 1
#     cols = len(seq2) + 1
#     score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
# # Initialize the traceback matrix
#     traceback_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

# # Initialize the first row and column of the scoring matrix
#     for i in range(1, rows):
#         score_matrix[i][0] = score_matrix[i-1][0] + gap
#         traceback_matrix[i][0] = 1

#     for j in range(1, cols):
#         score_matrix[0][j] = score_matrix[0][j-1] + gap
#         traceback_matrix[0][j] = 2

# # Fill in the scoring matrix
#     for i in range(1, rows):
#         for j in range(1, cols):
#             match_mismatch = match if seq1[i-1] == seq2[j-1] else mismatch
#             diag_score = score_matrix[i-1][j-1] + match_mismatch
#             up_score = score_matrix[i-1][j] + gap
#             left_score = score_matrix[i][j-1] + gap
#             max_score = max(diag_score, up_score, left_score)
#             score_matrix[i][j] = max_score

#     if max_score == diag_score:
#         traceback_matrix[i][j] = 3
#     elif max_score == up_score:
#         traceback_matrix[i][j] = 1
#     else:
#         traceback_matrix[i][j] = 2

# # Traceback to find the alignment
#     align1 = ""
#     align2 = ""
#     i, j = rows - 1, cols - 1
#     while i > 0 or j > 0:
#         if traceback_matrix[i][j] == 3:
#             align1 = seq1[i-1] + align1
#             align2 = seq2[j-1] + align2
#             i -= 1
#             j -= 1
#         elif traceback_matrix[i][j] == 1:
#             align1 = seq1[i-1] + align1
#             align2 = "-" + align2
#             i -= 1
#         else:
#             align1 = "-" + align1
#             align2 = seq2[j-1] + align2
#             j -= 1

#     return align1, align2

# # Example usage
# seq1 = "TACGCA"
# seq2 = "TATACGCA"
# alignment1, alignment2 = needleman_wunsch(seq1, seq2)
# print("Sequence 1:", alignment1)
# print("Sequence 2:", alignment2)
 


# 多序列比对

# def clustalw(sequences,match=1,mismatch=-1,gap=-1):
#     aligment = [list(seq)for seq in sequences]

#     # perform pairwise aligments
#     while len(aligment) >1:
#         scores = []
#         for i in range(len(aligment)):
#             for j in range(i+1,len(aligment)):
#                 score = sum(a == b for a,b in zip(aligment[i],aligment[j]))
#                 scores.append((i,j,score))
#                 i,j,_ = max(scores, key = lambda x:x[2])
#                 merged_seq=[]
#                 for a,b in zip(aligment[i],aligment[j]):
#                     if a==b:
#                         merged_seq.append(a)
#                     else:
#                         merged_seq.append("-")
#                         aligment[i] = merged_seq
#                         aligment.pop(j)
#                 return  ''.join(aligment[0])
# sequences = [
# "AGCT",
# "AGCT",
# "AGCT",
# "AGCT"
# ]
# alignment = clustalw(sequences)
# print("Multiple Sequence Alignment:")
# print(alignment)



# 生成随机序列
import random
def generate_random_sequence(length):
    bases = ["A","T","G","C"]
    sequence = ""
    for i in range(length):
        sequence += random.choice(bases)
    return sequence
sequence = generate_random_sequence(100)
# print(sequence)


# sortDNAsequencesbylength
def sortDNAsequencesbylength(sequences):
    sequences.sort(key = len)
    return sequences

# 比对多个序列
def clustalw(sequences, match=1, mismatch=-1, gap=-1):
    alignment = [list(seq) for seq in sequences]
    # Perform pairwise alignments
    while len(alignment) > 1:
        scores = []
        for i in range(len(alignment)):
            for j in range(i+1, len(alignment)):
                score = sum(a == b for a, b in zip(alignment[i], alignment[j]))
                scores.append((i, j, score))
        i, j, _ = max(scores, key=lambda x: x[2])
        merged_seq = []
        for a, b in zip(alignment[i], alignment[j]):
            if a == b:
                merged_seq.append(a)
            else:
                merged_seq.append("-")
        alignment[i] = merged_seq
        alignment.pop(j)
    return "".join(alignment[0])
# 使用generate_random_sequence函数生成50个随机序列
DNAsequence = [generate_random_sequence(random.randint(50,100)) for _ in range(50)]
# print(DNAsequence)
# print(sortDNAsequencesbylength(DNAsequence))
print(clustalw(DNAsequence))












