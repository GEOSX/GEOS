#!/usr/bin/env python
# Python wrapper script for generating the correct cmake line with the options specified by the user.
#
# Please keep parser option names as close to possible as the names of the cmake options they are wrapping.

import argparse
import logging
import os
import shutil
import stat
import subprocess
import sys


def extract_cmake_location(file_path):
    logging.info("Extracting cmake entry from host config file " + file_path)
    if os.path.exists(file_path):
        cmake_line_prefix = "# cmake executable path: "
        file_handle = open(file_path, "r")
        content = file_handle.readlines()
        for line in content:
            if line.startswith(cmake_line_prefix):
                return line.split(" ")[4].strip()
        logging.info("Could not find a cmake entry in host config file. Using ${PATH}.")
    return None


def parse_args(cli_arguments):
    """
    Parse command line arguments into an ArgumentParser instance.
    :param cli_arguments: The command line arguments as an array of string.
    :return: The ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(description="Configure cmake build.")

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-bp",
        "--buildpath",
        dest="build_path",
        type=str,
        default="",
        help=
        "Specify path for build directory. If both `--buildpath` and `--buildrootdir` are not specified, the build directory will be created in current directory, in a subdirectory named `build-<config-file-name>-<buildtype>`. The `--buildpath` option is not compatible with `--buildrootdir`."
    )

    group.add_argument(
        "-br",
        "--buildrootdir",
        dest="build_root_dir",
        type=str,
        default="",
        metavar="BUILD_ROOT_DIR",
        help=
        "Specify path for root of build directory. The build directory will be created as `BUILD_ROOT_DIR/build-<config-file-name>-<buildtype>`. The `--buildrootdir` option is not compatible with `--buildpath`."
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-ip",
        "--installpath",
        dest="install_path",
        type=str,
        default="",
        help=
        "Specify path for installation directory. If both `--installpath` and `--buildpathdir` are not specified, the install directory will be created in current directory, in a subdirectory named `install-<config-file-name>-<buildtype>`. The `--installpath` option is not compatible with `--installrootdir`."
    )

    group.add_argument(
        "-ir",
        "--installrootdir",
        dest="install_root_dir",
        type=str,
        default="",
        metavar="INSTALL_ROOT_DIR",
        help=
        "Specify path for root of install directory. The install directory will be created as `INSTALL_ROOT_DIR/build-<config-file-name>-<buildtype>`. The `--installrootdir` option is not compatible with `--installpath`."
    )

    parser.add_argument("--noinstall",
                        dest="no_install",
                        action="store_true",
                        help="Do not create an install directory.")

    parser.add_argument("-bt",
                        "--buildtype",
                        dest="build_type",
                        type=str,
                        choices=["Release", "Debug", "RelWithDebInfo", "MinSizeRel"],
                        default="Debug",
                        help="build type.")

    parser.add_argument("-e", "--eclipse", action='store_true', help="Create an eclipse project file.")
    parser.add_argument("-n", "--ninja", action='store_true', help="Create a ninja project.")
    parser.add_argument("-x", "--xcode", action='store_true', help="Create an xcode project.")

    parser.add_argument(
        "-ecc",
        "--exportcompilercommands",
        dest="export_compiler_commands",
        action='store_true',
        help=
        "generate a compilation database.  Can be used by the clang tools such as clang-modernize.  Will create a file called 'compile_commands.json' in build directory."
    )

    parser.add_argument("-hc",
                        "--hostconfig",
                        dest="host_config",
                        required=True,
                        type=str,
                        help="select a specific host-config file to initalize CMake's cache")

    parser.add_argument("-gvz", "--graphviz", action="store_true", help="Generate graphviz dependency graph")

    args, unknown_args = parser.parse_known_args(cli_arguments)
    if unknown_args:
        logging.info("Passing the following unknown arguments directly to cmake: %s" % unknown_args)
    return args, unknown_args


def main(calling_script, args, unknown_args):
    ########################
    # Find CMake Cache File
    ########################
    scripts_dir = os.path.dirname(os.path.abspath(calling_script))

    cache_file = os.path.abspath(args.host_config)
    platform_info = os.path.split(cache_file)[1]
    if platform_info.endswith(".cmake"):
        platform_info = platform_info[:-6]

    assert os.path.exists(cache_file), "Could not find cmake cache file '%s'." % cache_file
    logging.info("Using host config file: '%s'." % cache_file)

    #####################
    # Setup Build Dir
    #####################
    if args.build_path:
        # use explicit build path
        build_path = args.build_path
    else:
        # use platform info & build type
        build_path = "-".join(["build", platform_info, args.build_type.lower()])
        if args.build_root_dir != "":
            build_path = os.path.join(args.build_root_dir, build_path)

    logging.info("Build path is: " + build_path)

    build_path = os.path.abspath(build_path)

    if os.path.exists(build_path):
        logging.info("Build directory '%s' already exists. Deleting..." % build_path)
        shutil.rmtree(build_path)

    logging.info("Creating build directory '%s'..." % build_path)
    os.makedirs(build_path)

    #####################
    # Setup Install Dir
    #####################
    # For install directory, we will clean up old ones, but we don't need to create it, cmake will do that.
    if not args.no_install:
        if args.install_path != "":
            install_path = os.path.abspath(args.install_path)
        else:
            # use platform info & build type
            install_path = "-".join(["install", platform_info, args.build_type.lower()])
            if args.install_root_dir != "":
                install_path = os.path.join(args.install_root_dir, install_path)

        install_path = os.path.abspath(install_path)

        #     logging.info("Install directory '%s' already exists. Deleting..." % install_path)
        #     shutil.rmtree(install_path)

        if not os.path.exists(install_path):
            logging.info("Creating install path '%s'..." % install_path)
            os.makedirs(install_path)

    ############################
    # Build CMake command line
    ############################

    cmake_line = list()
    cmake_line.append(extract_cmake_location(cache_file) or "cmake")

    # Add build type (opt or debug)
    cmake_line.append("-DCMAKE_BUILD_TYPE=" + args.build_type)

    if not args.no_install:
        cmake_line.append("-DCMAKE_INSTALL_PREFIX=%s" % install_path)

    if args.export_compiler_commands:
        cmake_line.append("-DCMAKE_EXPORT_COMPILE_COMMANDS=ON")

    if args.eclipse:
        cmake_line.append('-G"Eclipse CDT4 - Unix Makefiles"')

    if args.ninja:
        cmake_line.append('-GNinja')

    if args.xcode:
        cmake_line.append('-GXcode')

    if args.graphviz:
        cmake_line.append("--graphviz=dependency.dot")
        dot_line = "dot -Tpng dependency.dot -o dependency.png"

    for unknown_arg in unknown_args:
        if not unknown_arg.startswith('-D'):
            logging.warning("Additional argument '%s' does not start with '-D'. Keeping it nevertheless." % unknown_arg)
        cmake_line.append(unknown_arg)

    # Append cache file at the end of the command line to make previous argument visible to the cache.
    cmake_line.append("-C%s" % cache_file)

    cmake_line.append(os.path.normpath(os.path.join(scripts_dir, "..", "src")))

    # Dump the cmake command to file for convenience
    cmake_cmd = os.path.join(build_path, "cmake_cmd")
    with open(cmake_cmd, "w") as cmd_file:
        cmd_file.write(" ".join(cmake_line) + os.linesep)
    st = os.stat(cmake_cmd)
    os.chmod(cmake_cmd, st.st_mode | stat.S_IEXEC)

    ############################
    # Run CMake
    ############################
    logging.info("Changing to build directory '%s'" % build_path)
    os.chdir(build_path)
    with open(cmake_cmd, "r") as cmd_file:
        logging.info("Executing cmake line: '%s'" % cmd_file.read().rstrip(os.linesep))
    subprocess.call(cmake_cmd, shell=True)

    if args.graphviz:
        subprocess.call(dot_line, shell=True)


if __name__ == '__main__':
    logging.basicConfig(format='[%(filename)s]:[%(levelname)s]: %(message)s', level=logging.INFO)
    main(sys.argv[0], *parse_args(sys.argv[1:]))
