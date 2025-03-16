#!/bin/bash

# 讀取 fastq_files.txt 中的每一行
while read FILE; do
  bgzip "$FILE"
done < fastq_files.txt
