#!/usr/bin/env python

###############################################################################
# Copyright (c) 2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-746361
#
# All rights reserved. See COPYRIGHT for details.
#
# This file is part of the GEOSX Simulation Framework.
#
# GEOSX is a free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the FREE
# Software Foundation) version 2.1 dated February 1999.
###############################################################################

# Python script to add GEOSX LLNL copyright notice at top of files in a directory
#
# The script takes a directory and checks all files in the directory.
# If the second line in the file does not match the copyright string, we add it.
#
# Modified from an initial script by P. Sinha

# this is a handy command to check if there are any files that have not be changed from HEAD
# git ls-files --full-name | grep -v "$(git diff --name-only HEAD)"

import os
import sys
import argparse


copyright_str = \
"""/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron 
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */
"""

old_copyright_str = \
"""
/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */
"""

copyright_str_arr = copyright_str.split("\n")[:-1]

old_copyright_str_arr = old_copyright_str.split("\n")[:-1]

max_copyright_lines = max(len(copyright_str_arr), len(old_copyright_str_arr))


def getLineToBeginAt(lines):
    """Given the LINES of a file return the line at which to begin writing the new copyright
       header. If no header needs to be written return -1.
    """

    # Check if the file has the copyright header.
    copyright_match = False
    if len(lines) >= len(copyright_str_arr):
        copyright_match = True

        for i in range(len(copyright_str_arr)):
            if i == 1:
                continue
            if lines[i][:-1] != copyright_str_arr[i]:
                copyright_match = False
                break

    # If the file has the copyright header check if it needs to be updated.
    if copyright_match:
        if lines[1][:-1] == copyright_str_arr[1]:
            print( "\t Already has copyright statement." )
            return -1

        print( "\t Need to update copyright statement." )
        return len(copyright_str_arr)

    # Check if the file has the old header.
    old_copyright_match = False
    if len(lines) >= len(old_copyright_str_arr):
        old_copyright_match = True

        for i in range(len(old_copyright_str_arr)):
            if lines[i][:-1] != old_copyright_str_arr[i]:
                old_copyright_match = False
                break

    if old_copyright_match:
        print( "\t Has old copyright statement." )
        return len(old_copyright_str_arr)

    print( "\t Missing copyright statement." )
    return -2


def checkAndAddCopyrightHeader(filename, testOnly=False):
    """Update the copyright header on the given file. If testOnly file is not updated.
    """

    with open(filename, "r+") as f:
        # print( "  Processing file {}:".format(filename) )

        lines = []
        for i in range(max_copyright_lines):
            line = f.readline()
            if line == "":
                break

            lines += [line]

        line_to_begin_at = getLineToBeginAt(lines)

        if line_to_begin_at == -2:
            print("Check missing or malformed copyright (including empty lines) : {}".format(filename))

        if line_to_begin_at >= 0 and not testOnly:
            lines += f.readlines()
            f.seek(0)
            f.write(copyright_str)
            f.writelines(lines[line_to_begin_at:])
            print ( "\t Prepended copyright statement." )


def fileNameGenerator(rootDir, validExtensions, isRecursive=False):
    """Generator function for file names whose extensions are in the validExtensions tuple.
       Files are rooted in rootDir, and process is recursive if isRecursive==True.
    """

    if isRecursive:
        for path, dirlist, filelist in os.walk(rootDir):
            for f in (f for f in filelist if f.lower().endswith(validExtensions)):
                yield os.path.join(path, f)
    else:
        for f in os.listdir(rootDir):
            if f.lower().endswith(validExtensions):
                yield os.path.join(rootDir, f)


if __name__ == "__main__":

    ## Setup the argument parser, dir is required first argument
    parser = argparse.ArgumentParser(description="Append LLNL copyright message to files.")
    parser.add_argument("dir", type=str, help="specify directory containing files on which we want to operate."
                        )    # TODO -- should we accept multiple directories?

    # Option to recursively search for files
    parser.add_argument("-r",
                        "--recursive",
                        dest='isRecursive',
                        action='store_true',
                        help="add flag to recursively descend to subdirectories.")
    parser.add_argument("--no-recursive", dest='isRecursive', action='store_false')
    parser.set_defaults(isRecursive=False)

    # Test run to see which files might require a copyright notice.
    parser.add_argument(
        "-t",
        "--test",
        action='store_true',
        help="add flag if we only want to see which files are missing copyright statements, but not to modify any files."
    )

    # Additional possible featues to be implemented
    #
    # Specify a collection of file types.
    #parser.add_argument("-f", "--filetypes", type=str, nargs="*", help="add file types on which to operate ")    # default should be *.hpp and *.cpp
    # For now, we hardcode to c/c++ header and source files
    valid_extensions = (".hpp", ".cpp", ".h", ".c", ".cxx", ".cc")

    args = parser.parse_args()

    ## Iterate through files, check for and add copyright notice
    print( "Looking at directory {}".format(args.dir) )
    for fullFileName in fileNameGenerator(args.dir, valid_extensions, args.isRecursive):
        checkAndAddCopyrightHeader(fullFileName, args.test)
