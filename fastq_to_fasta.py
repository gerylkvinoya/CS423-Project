# !/usr/bin/env python3
 
# -*- coding: utf-8 -*-
 
from pysam import FastxFile
import textwrap
import os

#global variable to hold the list of the sequences tuple
sequenceList = []

########################################################
# read_fasta_q_file -- parse fastq file using pysam.FastxFile
#
#   fasta_q_file - name of fastq file to read
#
#   return: the sequence id and sequence as a tuple
#   
#   source:
#       https://onestopdataanalysis.com/fastq-to-fasta/ 
#######################################################
def read_fasta_q_file(fasta_q_file):
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

    # the commented code creates a multi-fasta file
    # for s in sequenceList:
    #     f.write("> " + str(s[0]) + "\n")
    #     f.write(textwrap.fill(s[1], width=60))
    #     f.write("\n\n")

    combinedSequence = ""

    f.write("> combined sequence\n")
    for s in sequenceList:
        combinedSequence += trimSequence(s[1], 0, len(s))

    f.write(textwrap.fill(combinedSequence, width=60))


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

########################################################
# trimSequence -- RECURISVE
#                   trim the sequence based on the following primers
#                   sequences may have one or more of these:
#                   GCGGTAATTCCAGCTCCAATAG 
#                   CTCTGACAATGGAATACGAATA
#                   AAGGAGAAATHAATGTCT 
#                   AARCAACCTTGTGTAAGTCTC 
#
#   sequence - sequence to trim
#
#   return: trimmed sequence
#######################################################
def trimSequence(sequence, index, end):
    #BASE CASE: if we have reached the end of the sequence, return the sequence
    if index == end:
        return sequence

    newSeq = sequence

    #check for the two sequences that have length 22
    if len(sequence[index:len(sequence)]) >= 22:
        #check for sequence GCGGTAATTCCAGCTCCAATAG
        if sequence[index:index + 22] == "GCGGTAATTCCAGCTCCAATAG":
            newSeq = sequence[0:index] + sequence [index + 22:end]
            return trimSequence(newSeq, index, end)

        #check for sequence CTCTGACAATGGAATACGAATA
        if sequence[index:index + 22] == "CTCTGACAATGGAATACGAATA":
            newSeq = sequence[0:index] + sequence [index + 22:end]
            return trimSequence(newSeq, index, end) 

    #check for one sequence that has length 18 AAGGAGAAATHAATGTCT
    if len(sequence[index:len(sequence)]) >= 18:
        if sequence[index:index + 18] == "AAGGAGAAATHAATGTCT":
            newSeq = sequence[0:index] + sequence [index + 18:end]
            return trimSequence(newSeq, index, end)
    
    #check for one sequence that has length 21 AARCAACCTTGTGTAAGTCTC
    if len(sequence[index:len(sequence)]) >= 21:
        if sequence[index:index + 21] == "AARCAACCTTGTGTAAGTCTC":
            newSeq = sequence[0:index] + sequence [index + 21:end]
            return trimSequence(newSeq, index, end)
    
    index += 1
    return trimSequence(newSeq, index, end)

searchDir('sequences')
writeFasta("out.fasta", sequenceList)

#s = "QAARCAACCTTGTGTAAGTCTCAAGGAGAAATHAATGTCTQAAGGAGAAATHAATGTCTQCTCTGACAATGGAATACGAATACTCTGACAATGGAATACGAATAQGCGGTAATTCCAGCTCCAATAGQGCGGTAATTCCAGCTCCAATAGQGCGGTAATTCCAGCTCCAATAGAARCAACCTTGTGTAAGTCTCAAGGAGAAATHAATGTCT"

#print(trimSequence(s, 0, len(s)))
