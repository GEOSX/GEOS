#!/usr/bin/env python3
# Python script to
import sys
import os
import stat
import subprocess
import argparse
import platform
import shutil


# fines all files recursively from
def findFiles(folder, extension):
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if (extension in filename):
                yield os.path.join(root, filename)


parser = argparse.ArgumentParser(description='Process ats output to filter diffs.')

parser.add_argument('-d',
                    '--directory',
                    type=str,
                    default='integratedTests',
                    help='directory to search recursively for files with specified extension')

parser.add_argument('-ext', '--extension', type=str, default='.log', help='extension of files to filter')

parser.add_argument('-tl',
                    '--numTrailingLines',
                    type=int,
                    default=5,
                    help='number of lines to include in block after match is found.')

parser.add_argument('-e',
                    '--exclusionStrings',
                    type=str,
                    nargs="*",
                    default=[],
                    help='What stings to look for in order to exclude a block')

args, unknown_args = parser.parse_known_args()
if unknown_args:
    print("unknown arguments %s" % unknown_args)

# What strings to look for in order to flag a line/block for output
matchStrings = ['Error:']

# What stings to look for in order to exclude a block
#exclusionStrings = [ 'sizedFromParent', 'different shapes' ]
#exclusionStrings = [ 'sizedFromParent', 'different shapes', 'but not the' ]
exclusionStrings = ['logLevel', 'NonlinearSolverParameters', 'has a child', 'different shapes', 'different types', 'differing types']
exclusionStrings += args.exclusionStrings

directory = args.directory
extension = args.extension
numTrailingLines = args.numTrailingLines

for fileName in findFiles(directory, extension):

    filteredErrors = ''

    with open(fileName) as f:
        lines = f.readlines()

        for i in range(0, len(lines)):
            line = lines[i]
            if all(matchString in line for matchString in matchStrings):
                matchBlock = '  ' + lines[i - 1]
                matchBlock += '  ' + line

                for j in range(1, numTrailingLines + 1):
                    if i + j >= len(lines):
                        matchBlock += '  ***** No closing line. file truncated? Filters may not be properly applied! *****'
                        break
                    matchBlock += '  ' + lines[i + j]

                    if ('******************************************************************************'
                            in lines[i + j]):
                        break
                i += j

                if not any(excludeString in matchBlock for excludeString in exclusionStrings):
                    filteredErrors += matchBlock

    if (len(filteredErrors)):
        print("\nFound unfiltered diff in: ", fileName)
        print(filteredErrors, flush=True)

#for i in range(1,1+1):
#    print( i )
