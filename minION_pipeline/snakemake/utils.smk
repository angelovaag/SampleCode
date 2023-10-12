import os
import pandas as pd
import csv

def read_sampList_file():
    with open('sampList.txt') as f:
    	samples = f.read().splitlines()
    return samples
samples=read_sampList_file()


#samples=["ERR3201932", "ERR3201933","Mock", "Tmu_0"]

