#教程来自 https://biopython.org/docs/latest/index.html
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, SimpleLocation
import gzip

# My_seq = Seq("ATCCCGCTCGAGGAGGAGAGTTTCAGAGATCACGAATACATCCATATTACCCAGAGAGAG") # Create a sequence object
# # print(My_seq) # Print the sequence
# My_complement_seq = My_seq.complement() # Return the complement of the sequence



# for seq_record in SeqIO.parse("./biopython/Doc/examples/ls_orchid.fasta", "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# for seq_record in SeqIO.parse("./biopython/Doc/examples/ls_orchid.gbk", "genbank"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# for index, letter in enumerate(My_seq):
#     print("%i %s" % (index, letter))
# print(My_seq[25:50]) 

# n = "AAAA".count("AA")
# m =Seq("AAAA").count("AA")

# My_Seq_G_number = My_seq.count("G")
# My_Seq_C_number = My_seq.count("C")
# My_Seq_T_number = My_seq.count("T")
# My_Seq_A_number = My_seq.count("A")
# My_seq_CG_contend = 100 * (My_seq.count("G") + My_seq.count("C")) / len(My_seq)
# print(My_Seq_G_number)
# print(My_Seq_C_number)
# print(My_Seq_T_number)
# print(My_Seq_A_number)
# print(f"CG所占百分比:{My_seq_CG_contend:.2f}%")

# My_aeq_gc = gc_fraction(My_seq)
# print(f"GC所占百分比:{My_aeq_gc:.2f}")

# Slicing a sequence

# n = My_seq[4:12]
# m = My_seq[::-1]
# print(n)
# print(m)

# Turning Seq objects into strings
# n = str(My_seq)
# print(type(My_seq))
# print(type(n))
# fasta_format_string = ">Name seq1 \n%s\n" % My_seq
# print(fasta_format_string)

# Concatenating or adding sequences
# seq1 = Seq("ACGT")
# seq2 = Seq("AACCGG")
# seq3 = seq1 + seq2
# print(seq3)

# protein_seq = Seq("EVRNAK")
# dna_seq = Seq("ACGT")
# confusing_seq = protein_seq + dna_seq
# print(confusing_seq)

# list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
# concatenated = Seq("_")
# for s in list_of_seqs:
#     concatenated += s
# print(concatenated)

# contigs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
# spacer = Seq("N" * 10)
# print(spacer.join(contigs))


# Changing case
# dna_seq = Seq("acgtACGT")
# print(dna_seq)
# seq_up = dna_seq.upper()
# seq_low = dna_seq.lower()
# print(seq_up)
# print(seq_low)
# print("GTAC" in dna_seq)

# Nucleotide sequences and (reverse) complements
# my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
# print(my_seq)
# my_seq_com = my_seq.complement()
# my_seq_recom = my_seq.reverse_complement()
# print(my_seq_com)
# print(my_seq_recom)
# print(my_seq[::-1]) #与原数据顺序相反


# Transcription
# coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
# print(coding_dna)
# template_dna = coding_dna.reverse_complement()
# print(template_dna)
# messenger_rna = coding_dna.transcribe()
# print(messenger_rna)
# print(template_dna.reverse_complement().transcribe())

# messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
# print(Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG'))
# print(messenger_rna.back_transcribe())


# Translation
# messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
# print(messenger_rna)
# protein = messenger_rna.translate()
# print(protein)

# coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
# print(coding_dna)
# coding_dna.translate()
# print(coding_dna.translate(table="Vertebrate Mitochondrial"))
# print(coding_dna.translate(table=2))

# print(coding_dna.translate(),
# coding_dna.translate(to_stop=True),
# coding_dna.translate(table=2),
# coding_dna.translate(table=2, to_stop=True),
# coding_dna.translate(table=2, stop_symbol="@"))


#例子
# gene = Seq(
#     "GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA"
#     "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT"
#     "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT"
#     "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT"
#     "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA"
# )
# protein =gene.translate(table="Bacterial",to_stop=True)
# protein2 =gene.translate(table="Bacterial", cds=True)
# print(protein2,protein)


# Translation Tables
# standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
# print(standard_table)
# mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
# print(mito_table)
# standard_table = CodonTable.unambiguous_dna_by_id[1]
# mito_table = CodonTable.unambiguous_dna_by_id[2]

# stop_condons = mito_table.stop_codons
# print(stop_condons)
# start_condons = mito_table.start_codons
# ch = mito_table.forward_table["ACG"]


# # Comparing Seq objects
# seq1 = Seq("ACGT")
# print("ACGT" == seq1) #True


# # Sequences with unknown sequence contents
# unknown_seq = Seq(None, 10)
# len(unknown_seq)


# Sequences with partially defined sequence contents
# seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length=159345973)
# print(seq[117512685:117512690])

# seq = Seq("ACGT")
# undefined_seq = Seq("AAAAAAAAAAAA", length=10)
# print(seq + undefined_seq + seq)


# MutableSeq objects
# my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
# # my_seq[5] = "G"
# # my_seq = my_seq[:5] + "G" + my_seq[6:] #正确修改方式
# mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
# # print(my_seq)
# # print(mutable_seq)
# mutable_seq[5] = "C"
# mutable_seq.remove("T")
# mutable_seq.reverse()
# new_seq = Seq(mutable_seq) #将编译序列改编成Seq格式


# # Finding subsequences
# seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
# seq.index("ATGGGCCGC")
# seq.index(b"ATGGGCCGC")
# seq.index(bytearray(b"ATGGGCCGC"))
# seq.index(Seq("ATGGGCCGC"))
# seq.index(MutableSeq("ATGGGCCGC"))
# # seq.index("ACTG")  #没有的话就会报错
# # seq.find("ACTG") #-1

# for index, sub in seq.search(["CC", "GGG", "CC"]):
#     print(index, sub)


# # Working with strings directly
# from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
# my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
# reverse_complement(my_string)
# transcribe(my_string)
# back_transcribe(my_string)
# translate(my_string)






# Sequence annotation objects
## 使用的是SeqRecord类
### 如何导入SeqRecord类
# # from Bio import SeqIO
# for index, record in enumerate(SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank")):
#     print(
#         "index %i, ID = %s, length %i, with %i features"
#         % (index, record.id, len(record.seq), len(record.features))
#     )
# print(record)
# print(dir(record))
# print(record.seq)
# print(record.annotations)
# print(record.format("fasta"))

### 创建一个SeqRecord对象
# simple_seq = Seq("GATC")
# simple_seq_r = SeqRecord(simple_seq)
# simple_seq_r.name= "SimpleSeq"
# simple_seq_r.id = "AC12345"
# simple_seq_r.description = "Made up sequence I wish I could write a paper about"
# simple_seq_r.annotations["evidence"] = " None. I just made it up."
# simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
# simple_seq_r = SeqRecord(simple_seq, 
#                          id="AC12345",
#                          name = "SimpleSeq", 
#                          description="Made up sequence I wish I could write a paper about",
#                          annotations= {"evidence" : " None. I just made it up."},
#                          letter_annotations= {"phred_quality" : [40, 40, 38, 30]})
# print(simple_seq_r)
# print(simple_seq)
# print(simple_seq_r.annotations["evidence"])
# print(simple_seq_r.letter_annotations["phred_quality"])


# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

# record = SeqRecord(
#     Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
#     id="YP_025292.1",
#     name="HokC",
#     description="toxic membrane protein, small",
# )
# print(record.name)
# print(record.description)
# SeqIO.write(record, "test.fasta", "fasta")

# # SeqRecord objects from FASTA files
# record = SeqIO.read("NC_005816.fna", "fasta")
# record_1 = SeqIO.read("NC_005816.gb", "genbank")
# # print(record.id)
# My_seq = Seq(record.seq)
# # print(My_seq)
# print(type(record))
# print(type(My_seq))


## Feature, location and position objects
# start_pos = SeqFeature.AfterPosition(5)
# end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
# my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
# print(my_location)
# print(int(my_location.start))
# print(int(my_location.end))
# exact_location = SeqFeature.SimpleLocation(5, 9)
# print(exact_location)

# my_snp = 4350
# record = SeqIO.read("NC_005816.gb", "genbank")
# for feature in record.features:
#     if my_snp in feature:
#         print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))

# seq = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
# feature = SeqFeature(SimpleLocation(5, 18, strand=-1), type="gene")
# feature_seq = seq[feature.location.start : feature.location.end].reverse_complement()
# print(feature_seq)
# feature_seq = feature.extract(seq)
# print(feature_seq)



## comparsion
# record1 = SeqRecord(Seq("ACGT"), id="test")
# record2 = SeqRecord(Seq("ACGT"), id="test")
# n = record1.id == record2.id
# m = record2.seq == record1.seq
# print(n == m)


# # The format method
# record = SeqRecord(
#     Seq(
#         "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
#         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
#         "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
#         "SSAC"
#     ),
#     id="gi|14150838|gb|AAK54648.1|AF376133_1",
#     description="chalcone synthase [Cucumis sativus]",
# )
# print(record.format("fasta"))


# # Slicing a SeqRecord
# record = SeqIO.read("NC_005816.gb", "genbank")
# record
# len(record)
# len(record.features)
# print(record.features[20])
# sub_record = record[4300:4800]



# Adding SeqRecord objects
# record = next(SeqIO.parse("example.fastq", "fastq"))
# rec = SeqIO.read("NC_005816.gb", "genbank")
# print(rec.id, len(rec), len(rec.features), len(rec.dbxrefs), len(rec.annotations))





# Sequence Input/Output
# from Bio import SeqIO
# help(SeqIO)
## Parsing or Reading Sequences
## Reading Sequence Files
# for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))


# for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))
                       #||
# identifiers = [seq_record.id for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank")]
# print(identifiers)

## Iterating over the records in a sequence file
# record_iterator = SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta")
# first_record = next(record_iterator)
# print(first_record.id)
# print(first_record.description)

# second_record = next(record_iterator)
# print(second_record.id)
# print(second_record.description)

# Third_record = next(record_iterator)
# print(Third_record.id)
# print(Third_record.description)

# first_record = next(SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank"))

## Getting a list of the records in a sequence file
# records = list(SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank"))
# print("Found %i records" % len(records))
# print("The last record")
# last_record = records[-1]  # using Python's list tricks
# print(last_record.id)
# print(repr(last_record.seq))
# print(len(last_record))
# print("The first record")
# first_record = records[0]  # remember, Python counts from zero
# print(first_record.id)
# print(repr(first_record.seq))
# print(len(first_record))


## Extracting data
# record_iterator = SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank")
# first_record = next(record_iterator)
# print(first_record.annotations.keys())
# print(first_record.annotations.values())
# print(first_record.annotations['source'])

# all_species = []
# for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank"):
#     all_species.append(seq_record.annotations["organism"])
# all_species = [
#     seq_record.annotations["organism"]
#     for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank")
# ]

# for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta"):
#     print(seq_record.description)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#     all_species.append(seq_record.description.split()[1])

# all_species =[
#     seq_record.description.split()[1]
#     for seq_record in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta")]
# print(all_species)


## Modifying data
# record_iterator = SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta")
# first_record = next(record_iterator)
# print(first_record.id)
# first_record.id = "new_id"
# print(first_record.id)
# first_record.description = first_record.id + " " + "desired new description"
# print(first_record.format("fasta")[:200])


## Parsing sequences from compressed files
# print(sum(len(r) for r in SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "gb")))
# with open("./Python_learning/biopython/Doc/examples/ls_orchid.gbk") as handle:
#     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))

# with gzip.open("./Python_learning/biopython/Doc/examples/ls_orchid.gbk.gz", "rt") as handle:
#     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))


## Parsing sequences from the net
#Note that just because you can download sequence data and parse it into a SeqRecord object in one go doesn’t mean this is a good idea. In general, you should probably download sequences once and save them to a file for reuse.

## Parsing GenBank records from the net
# from Bio import Entrez
# Entrez.email = "858180954@qq.com"
# with Entrez.efetch(
#     db="nucleotide", rettype="fasta", retmode="text", id="6273291"
# ) as handle:
#     seq_record = SeqIO.read(handle, "fasta")
# print("%s with %i features" % (seq_record.id, len(seq_record.features)))

# with Entrez.efetch(
#     db="nucleotide", rettype="gb", retmode="text", id="6273291"
# ) as handle:
#     seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"
# print("%s with %i features" % (seq_record.id, len(seq_record.features)))

# with Entrez.efetch(
#     db="nucleotide", rettype="gb", retmode="text", id="6273291,6273290,6273289"
# ) as handle:
#     for seq_record in SeqIO.parse(handle, "gb"):
#         print("%s %s..." % (seq_record.id, seq_record.description[:50]))
#         print(
#             "Sequence length %i, %i features, from: %s"
#             % (
#                 len(seq_record),
#                 len(seq_record.features),
#                 seq_record.annotations["source"],
#             )
#         )


## Parsing SwissProt sequences from the net
# from Bio import ExPASy
# with ExPASy.get_sprot_raw("O23729") as handle:
#     seq_record = SeqIO.read(handle, "swiss")
# print(seq_record.id)
# print(seq_record.name)
# print(seq_record.description)
# print(repr(seq_record.seq))
# print("Length %i" % len(seq_record))
# print(seq_record.annotations["keywords"])



## Sequence files as Dictionaries
# handle = open("./Python_learning/biopython/Tests/TwoBit/sequence.bigendian.2bit", "rb")
# records = SeqIO.parse(handle, "twobit")
# print(records.keys())
# print(records["seq222"].name)
# print(records["seq222"].seq)
# handle.close()


## Sequence files as Dictionaries – In memory
# orchid_dict = SeqIO.to_dict(SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank"))
# print(orchid_dict.keys())
# print(orchid_dict["Z78501.1"].seq)
# seq_record = orchid_dict["Z78475.1"]
# print(seq_record.description)


## Specifying the dictionary keys
# orchid_dict = SeqIO.to_dict(SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta"))
# print(orchid_dict.keys())
# def get_accession(record):
#     """Given a SeqRecord, return the accession number as a string.
#     e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
#     """
#     parts = record.id.split("|")
#     assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
#     return parts[3]
# orchid_dict = SeqIO.to_dict(
#     SeqIO.parse("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta"), key_function=get_accession
# )
# print(orchid_dict.keys())

## Indexing a dictionary using the SEGUID checksum
# from Bio.SeqUtils.CheckSum import seguid
# for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
#     print(record.id, seguid(record.seq))

# seguid_dict = SeqIO.to_dict(
#     SeqIO.parse("ls_orchid.gbk", "genbank"), lambda rec: seguid(rec.seq)
# )
# record = seguid_dict["MN/s0q9zDoCVEEc+k/IFwCNF2pY"]
# print(record.id)
# print(record.description)



## Sequence files as Dictionaries – Indexed files
# orchid_dict = SeqIO.index("ls_orchid.gbk", "genbank")
# print(len(orchid_dict))
# print(orchid_dict.keys())
# seq_record = orchid_dict["Z78475.1"]
# print(seq_record.description)
# orchid_dict.close()

# orchid_dict = SeqIO.index("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta")
# print(len(orchid_dict))
# print(orchid_dict.keys())
# def get_acc(identifier):
#     """Given a SeqRecord identifier string, return the accession number as a string.

#     e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
#     """
#     parts = identifier.split("|")
#     assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
#     return parts[3]
# orchid_dict = SeqIO.index("./Python_learning/biopython/Doc/examples/ls_orchid.fasta", "fasta", key_function=get_acc)
# print(orchid_dict.keys())



## Getting the raw data for a record
# uniprot = SeqIO.index("uniprot_sprot.dat", "swiss")
# with open("selected.dat", "wb") as out_handle:
#     for acc in ["P33487", "P19801", "P13689", "Q8JZQ5", "Q9TRC7"]:
#         out_handle.write(uniprot.get_raw(acc))



## Sequence files as Dictionaries – Database indexed files
# import glob
# from Bio import SeqIO
# files = glob.glob("gbvrl*.seq")
# print("%i files to index" % len(files))
# gb_vrl = SeqIO.index_db("gbvrl.idx", files, "genbank")
# print("%i sequences indexed" % len(gb_vrl))
# print(gb_vrl["AB811634.1"].description)
# print(gb_vrl.get_raw("AB811634.1"))




## Indexing compressed files
# orchid_dict = SeqIO.index("./Python_learning/biopython/Doc/examples/ls_orchid.gbk", "genbank")
# print(len(orchid_dict))
# print(orchid_dict.keys())
# orchid_dict = SeqIO.index("ls_orchid.gbk.bgz", "genbank")
# print(len(orchid_dict))





## Writing Sequence Files
# rec1 = SeqRecord(
#     Seq(
#         "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
#         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
#         "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
#         "SSAC",
#     ),
#     id="gi|14150838|gb|AAK54648.1|AF376133_1",
#     description="chalcone synthase [Cucumis sativus]",
# )

# rec2 = SeqRecord(
#     Seq(
#         "YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
#         "DMVVVEIPKLGKEAAVKAIKEWGQ",
#     ),
#     id="gi|13919613|gb|AAK33142.1|",
#     description="chalcone synthase [Fragaria vesca subsp. bracteata]",
# )

# rec3 = SeqRecord(
#     Seq(
#         "MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
#         "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
#         "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
#         "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
#         "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
#         "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
#         "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
#     ),
#     id="gi|13925890|gb|AAK49457.1|",
#     description="chalcone synthase [Nicotiana tabacum]",
# )
# my_records = [rec1, rec2, rec3]

# SeqIO.write(my_records, "my_example.faa", "fasta")


# records = SeqIO.parse("ls_orchid.gbk", "genbank")
# count = SeqIO.write(records, "my_example.fasta", "fasta")
# print("Converted %i records" % count)
# count = SeqIO.convert("ls_orchid.gbk", "genbank", "my_example.fasta", "fasta")
# print("Converted %i records" % count)

# for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
#     print(record.id)
#     print(record.seq.reverse_complement())

# records = [
#     rec.reverse_complement(id="rc_" + rec.id, description="reverse complement")
#     for rec in SeqIO.parse("ls_orchid.fasta", "fasta")
# ]

# records = [
#     rec.reverse_complement(id="rc_" + rec.id, description="reverse complement")
#     for rec in SeqIO.parse("ls_orchid.fasta", "fasta")
#     if len(rec) < 700
# ]
# SeqIO.write(records, "rev_comp.fasta", "fasta")


## Getting your SeqRecord objects as formatted strings
# from io import StringIO

# records = SeqIO.parse("ls_orchid.gbk", "genbank")
# out_handle = StringIO()
# SeqIO.write(records, out_handle, "fasta")
# fasta_data = out_handle.getvalue()
# print(fasta_data)


# with open("ls_orchid_long.tab", "w") as out_handle:
#     for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
#         if len(record) > 100:
#             out_handle.write(record.format("tab"))





## 2024/10/6
## Sequence alignments
import numpy as np
from Bio.Align import Alignment
# seqA = "CCGGTTTTT"
# seqB = "AGTTTAA"
# seqC = "AGGTTT"
# sequences = [seqA, seqB, seqC]
# coordinates = np.array([[1, 3, 4, 7, 9], [0, 2, 2, 5, 5], [0, 2, 3, 6, 6]])
# alignment = Alignment(sequences, coordinates)
# print(alignment.sequences)
# print(alignment.coordinates)


# lines = ["CGGTTTTT", "AG-TTT--", "AGGTTT--"]
# for line in lines:
#     print(line)
# lines = [line.encode() for line in lines]  # convert to bytes
# sequences, coordinates = Alignment.parse_printed_alignment(lines)
# sequences = [sequence.decode() for sequence in sequences]
# print(coordinates)
# sequences[0] = "C" + sequences[0]
# sequences[1] = sequences[1] + "AA"
# print(sequences)
# coordinates[0, :] += 1
# print(coordinates)
# alignment = Alignment(sequences, coordinates)
# print(alignment)

## 比对之后没有继续学习，只学习了下面的构建系统发育树的方法







