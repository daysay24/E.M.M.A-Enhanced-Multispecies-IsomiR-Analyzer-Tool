import pandas as pd 

seq = ""

with open('./bruteforcesearch.txt') as file:
    for line in file:
        seq += line.rstrip()
seq = seq.lower()
print(seq[27644:27675])