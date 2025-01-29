#!/bin/bash

# check if -a or --all is provided as 1st arg
METHOD=all
case $1 in
    -a|--all) METHOD=all; shift
    ;;
    -g|--git) METHOD=git; shift
    ;;
    *)
esac

# nothing to do if schema file not given
if [ -z "$1" ]; then
    echo "Usage: $0 [-a|--all|-g|--git] <schema> [<path>...]"
    exit
fi

SCHEMA=$1; shift
LOGFILE=xml_validation_results.log

# "-r" in GNU xargs omits the call if input is empty
# OS X xargs does not support it, but does the same by default
if [ "$(uname)" == "Darwin" ]; then
    XARGS="xargs"
else
    XARGS="xargs -r"
fi

# check if xmllint is present
if ! hash xmllint &> /dev/null; then
    >&2 echo "Error: xmllint is required to validate xml files"
    exit
fi

# check if git is present, if needed
if [ "$METHOD" = "git" ] && ! (hash git &> /dev/null); then
    >&2 echo "Error: git is required when -g or --git is specified"
    exit
else
    GIT_REPO_DIR=$1;
    echo "INFO: The specified GIT_REPO_DIR"
    echo $GIT_REPO_DIR
    # Check GIT_REPO_DIR and its existence. Then add it to Git safe directories
    if [ -z "$GIT_REPO_DIR" ]; then
        echo "Error: GIT_REPO_DIR is not set." && exit 1
    elif [ ! -d "$GIT_REPO_DIR" ]; then
        echo "Error: The specified GIT_REPO_DIR '$GIT_REPO_DIR' does not exist." && exit 1
    else
        echo "Adding '$GIT_REPO_DIR' to Git safe directories..."
        git config --global --add safe.directory "$GIT_REPO_DIR" && \
        echo "'$GIT_REPO_DIR' has been successfully added." || \
        echo "Error: Failed to add '$GIT_REPO_DIR'."
        echo "INFO: Printing safe.directory: "
        git config --global --list | grep safe.directory
    fi
fi

abs_path ()
{
    if [ "$#" -gt 0 ]; then
        realpath -s "$@"
    fi
}

list_xml_files_all () 
{
    abs_path $(find $1 -name "*.xml" -not -path "*/\.*")
}

# git does not have -C flag prior to 1.8.5, so we emulate it with cd
# this is actually the most reliable way, although not the fastest
# run in a subprocess so as to avoid having to cd back
list_xml_files_git ()
{
    local git_root=$(cd $path; git rev-parse --show-toplevel 2>/dev/null)
    if ! [ $? -eq 0 ]; then
        >&2 echo "Error: $path does not appear to be part of a git repository"
        exit 1
    fi
    local prefix=$(cd $path; git rev-parse --show-prefix 2>/dev/null)
    git --git-dir=$git_root"/.git" ls-files $prefix | grep -e ".*[.]xml$" | sed "s|^|$git_root/|g"
}


git --version
# emit location
pwd

# create/nullify the log file
echo -n > $LOGFILE

# validate each path separately and write results in the log
for path in "$@"; do
    # emit location
    echo $path
    cd "$path" || { echo "Directory not found: $path"; exit 1; }
    pwd
    ls -la
    git status
    git_root_c=$(git rev-parse --show-toplevel 2>/dev/null)
    prefix_c=$(git rev-parse --show-prefix 2>/dev/null)
    echo $git_root_c
    echo $prefix_c
    echo $git_root"/.git"
    collected_xml_files=$(list_xml_files_$METHOD $path)
    echo $collected_xml_files
    $collected_xml_files | $XARGS xmllint --schema $SCHEMA --noout >> $LOGFILE 2>&1
done

# print any failed validations on the stderr
grep -v validates $LOGFILE >&2

# if there are failed validations, message and return the status
if grep -q -v validates $LOGFILE; then
    >&2 echo "XML validation failed. See details in $LOGFILE"
    exit 1
fi

