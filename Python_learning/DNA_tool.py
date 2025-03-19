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
    sequence = input("请输入DNA序列:")
    analyze_tool = DNATool()
    analyze_tool.analyze(sequence)

if __name__ == "__main__":
    print("This will run only if the script is executed directly.")
    main()
