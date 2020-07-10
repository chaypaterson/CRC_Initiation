#!/bin/python
# Python 2.7 standard (but should be agnostic/extensible)
# This is a module that provides useful functions for analysing the data
# produced from the lineage tracking simulations


print("Importing pandas...")
import pandas
print("Done.\nLoading handler functions...")

# Selects all columns from the key file that match a given search key
# and outputs the column numbers

def KeyToCols(targetFile, Key):
    KeyFile = pandas.read_csv(targetFile,nrows=1)
    ColNums = []
    n = 0
    # Seek pattern "Key" in target key file
    for entry in KeyFile:
        # if entry contains 331, print column
        if Key in entry:
            ColNums.append(n)
        n+=1
    return ColNums

def StoreKeys(targetFile, Key):
    KeyFile = pandas.read_csv(targetFile,nrows=1)
    KeyList = []
    # Seek pattern "Key" in target key file
    for entry in KeyFile:
        # if entry contains 331, print column
        if Key in entry:
            KeyList.append(entry)
    return KeyList
    
# Selects matching columns from a target data file and prints them
def PrintCols(DataFile, ColSelect):
    TargetData = pandas.read_csv(DataFile)
    ColsFromFile=[]
    for i in ColSelect:
        ThisCol = TargetData.iloc[:,i]
        print(ThisCol.iloc[-1])
        ColsFromFile.append(ThisCol)

    return ColsFromFile

# Which of these paths are "hits"? And have non-zero probability?
def FindHits(DataFile, ColSelect):
    TargetData = pandas.read_csv(DataFile)
    HitCols=[]
    n = 0
    for i in ColSelect:
        ThisCol = TargetData.iloc[:,i]
        lastProb= ThisCol.iloc[-1]
        if (lastProb > 0):
            print(lastProb,n)
            HitCols.append(n)
        n += 1
    return HitCols

# Convert this list of column hits into a list of column keys ready to input
# into Mathematica (which numbers arrays from 1):
def MathematicColmHits(ColumHits, ColumKeys):
    MathematicaColmHits=[]
    for i in ColumHits:
        print(ColumKeys[i]+1)
        MathematicaColmHits.append(ColumKeys[i]+1)

    return MathematicaColmHits
