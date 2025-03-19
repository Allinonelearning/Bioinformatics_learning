from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os
from DNA_simple_analyse import DNAStructure

class DNALocator:
    """
    DNA定位
    """
    def __init__(self, dna_sequence):
        self.seq = Seq(dna_sequence)
    
    def find_motif(self, motif):
        positions = []
        start = 0
        while True:
            pos = self.seq.find(motif, start)
            
            if pos == -1:  # 如果motif未找到,退出循环
                break
            
            positions.append(pos)
            start = pos + 1
        
        return positions
    
    def analyze(self, motif):
        print(f"\n在DNA序列中查找 '{motif}' 的位置:")
        
        positions = self.find_motif(motif)
        
        if positions:
            print("发现以下位置:")
            for i, pos in enumerate(positions):
                print(f"{i+1}. 序列 '{motif}' 在位置 {pos+1}")
            
            # 保存结果到文件
            self.save_results(motif, positions)
        else:
            print("未发现指定的序列。")

    def save_results(self, motif, positions):
        file_name = f"{motif}_positions.txt"
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, file_name)

        with open(file_path, 'w') as f:
            f.write(f"在DNA序列中查找 '{motif}' 的位置:\n")
            f.write("发现以下位置:\n")
            for i, pos in enumerate(positions):
                f.write(f"{i+1}. 序列 '{motif}' 在位置 {pos+1}\n")

def main():
    dna_sequence = input("请输入DNA序列(只包含A、C、G和T): ").upper()
    
    while not set(dna_sequence).issubset('ACGT'):
        print("错误: 输入只包含A、C、G和T。请重新输入DNA序列:")
        dna_sequence = input().upper()
    
    locator = DNALocator(dna_sequence)
    
    while True:
        motif = input("\n请输入要查找的特定序列(或输入 'exit' 来退出): ").upper()
        
        if motif.lower() == 'exit':
            print("感谢使用 DNA 序列位置查找工具。")
            break
        
        locator.analyze(motif)

if __name__ == "__main__":
    main()

