from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

class DNAStructure:
    def __init__(self, sequence):
        self.seq = Seq(sequence)
    
    def reverse_complement(self):
        return self.seq.reverse_complement()
    
    def transcribe(self):
        return self.seq.transcribe()
    
    def translate(self, table="Standard"):
        return self.seq.translate(table=table)
    
    def gc_content(self):
        return gc_fraction(self.seq)

class DNATool:
    @staticmethod
    def analyze(sequence):
        dna = DNAStructure(sequence)
        print("DNA序列:", dna.seq)
        print("反向互补:", dna.reverse_complement())
        print("转录:", dna.transcribe())
        print("\n翻译结果:")
        print("标准密码子:", dna.translate())
        print("线粒体密码子:", dna.translate(table="Vertebrate Mitochondrial"))
        print("\nGC含量:", dna.gc_content())

def main():
    #sequence = input("请输入DNA序列:")
    while True:
        sequence = input("请输入DNA序列 (或输入 'exit' 来退出): ")
        
        if sequence.lower() == 'exit':
            print("感谢使用 DNA 序列统计工具。")
            break
        
        if not sequence or not all(base in 'ACGTacgt' for base in sequence):
            print("错误: 无效的输入。请确保输入只包含 A、C、G 和 T 字母（大小写不敏感）。")
            continue
        try:
            analyze_tool = DNATool()
            analyze_tool.analyze(sequence)
        except Exception as e:
            print(f"发生了一个未预期的错误: {str(e)}")

if __name__ == "__main__":
    main()






