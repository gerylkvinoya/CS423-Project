# !/usr/bin/env python3
 
# -*- coding: utf-8 -*-
 
from pysam import FastxFile
import textwrap
import os

#global variable to hold the list of the sequences tuple
sequenceList = []

########################################################
# read_fasta_q_file -- print as Fasta format
#
#   fasta_q_file - name of fastq file to read
#
#   return: the sequence id and sequence as a tuple
#   
#   source:
#       https://onestopdataanalysis.com/fastq-to-fasta/ 
#######################################################
def read_fasta_q_file(fasta_q_file):
    """Parse FASTA/Q file using `pysam.FastxFile`.
 
    Args:
 
        fasta_q_file (str): Path to FASTA/Q file.
 
    """
    with FastxFile(fasta_q_file) as fh:
        for entry in fh:
            sequence_id = entry.name
            sequence = entry.sequence
    
    return (sequence_id, sequence)

########################################################
# writeFasta -- write to file as Fasta format
#
#   filename - name of file to write to
#   seq_id - sequence id
#   seq - the sequence as a string
#######################################################
def writeFasta(filename, sequenceList):
    f = open(filename, "w")
    for s in sequenceList:
        f.write("> " + str(s[0]) + "\n")
        f.write(textwrap.fill(s[1], width=60))
        f.write("\n\n")
    f.close()

########################################################
# searchDir -- search the directory for fastq files
#
#   dirname - name of directory to search
#
#   return: list of sequences to write to the file
#######################################################
def searchDir(dirname):
    directory = dirname
    
    for file in os.listdir(directory):
        f = os.path.join(directory, file)

        #if directory, search that directory too
        if os.path.isdir(f):
            searchDir(f)

        if os.path.isfile(f):
            if f.endswith('.fastq'):
                sequenceList.append(read_fasta_q_file(f))


searchDir('sequences')
writeFasta("out.fasta", sequenceList)


