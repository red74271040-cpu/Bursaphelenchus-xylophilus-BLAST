import os
import json
import re
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def prepare_json_data(f_in):
    with open(f_in, "r") as f:
        return json.load(f)

def get_table_data(f_in):
    query = prepare_json_data(f_in)
    hit_counter = Counter(player['hit_name'] for player in query if player['hit_name'])
    efficient_counter = Counter(player['hit_name'] for player in query if player['is_efficient'] and player['hit_name'])
    
    table_data = []
    for hit_name, total_count in hit_counter.items():
        eff_count = efficient_counter.get(hit_name, 0)
        table_data.append([hit_name, total_count, eff_count])
    return table_data

def create_gbk(main_target_dict, off_target_dict, seq_file, out_file):
    # Genbank 생성 로직 (Python 3 대응)
    data = SeqIO.read(seq_file, "fasta")
    record = SeqRecord(data.seq, id=data.id, name=data.name, description=data.description)
    # 여기에 Feature 추가 로직 작성...
    with open(out_file, 'w') as f:
        SeqIO.write(record, f, "genbank")