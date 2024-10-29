#!/bin/bash

# Submodules not checking for
declare -ar exclusion_list=( "blt" "integratedTests" "uberenv" )
echo "Submodules that are excluded from sync test : ${exclusion_list[@]}"

# Do not pull large files
git lfs uninstall &> /dev/null

# Pull submodule to get .git files.
#git submodule update --init integratedTests
git submodule update --init src/cmake/blt
git submodule update --init src/coreComponents/LvArray
git submodule update --init src/coreComponents/constitutive/PVTPackage
git submodule update --init src/coreComponents/fileIO/coupling/hdf5_interface


# Initialize PR submodule hashes
declare -ar pr_hashes_array=( $(git submodule status | awk '{print $1}') )

# Initialize submodule paths
declare -ar paths_array=( $(git submodule status | awk '{print $2}') )

# Initialize differences between PR and origin/develop branches
declare -ar diff_array=( $(git diff --name-only origin/develop) )

# Initialize main branches for submodules
declare -Ar main_branches=(
  ["blt"]="origin/develop"
  ["LvArray"]="origin/develop"
  ["integratedTests"]="origin/develop"
  ["hdf5_interface"]="origin/master"
  ["PVTPackage"]="origin/master"
)


length=${#paths_array[@]}

# Returns exit code 0 if the hash for every submodule in the PR is equal
# to the hash of each submodule's main branch (develop or master).
# Note: For LvArray submodule only, checks if the PR hash is in the
#       main branch's commit history.
# Note: See "exclusion_list" for submodules that are exempted.
exit_code=0
unsync_submodules=()

for (( i=0; i<$length; i++))
do
  # Just the submodule name
  module_name="$(basename ${paths_array[$i]})"
  
  # Check if submodule is excluded from check
  excluded=0
  for ex in "${exclusion_list[@]}"
  do
    if [ "$module_name" = "$ex" ]
    then
      excluded=1
    fi
  done

  if [ $excluded -eq 1 ]
  then
    continue
  fi

  # Check hashes if not excluded

  # Submodule's main branch
  main_branch="${main_branches[$module_name]}"

  # Submodule main hash
  main_hash="$( git -C ${paths_array[$i]} rev-parse $main_branch )"

  # PR hash with prefixed character removed
  pr_hash="$( echo ${pr_hashes_array[$i]} | tr -cd [:alnum:] )"

  if [ $module_name == "LvArray" ]
  then
    lv_array_commits="$(git -C ${paths_array[$i]} log $main_branch --pretty=%H)"
    if grep -q "$pr_hash" <<< "$lv_array_commits"
    then
      echo "Submodule LvArray's main branch $main_branch contains"\
            "PR branch's latest commit : $pr_hash"
    else
      echo "Submodule LvArray's main branch $main_branch does not contain"\
            "PR branch's latest commit : $pr_hash"
      echo "---- PR branch has latest hash $pr_hash"
      echo "---- $module_name/$main_branch has latest hash $main_hash"
      unsync_submodules+=( "$module_name" )
      exit_code=1
    fi
  elif [ $pr_hash == $main_hash ]
  then
    echo "PR branch and $main_branch have the same hashes for submodule"\
         "$module_name : $pr_hash"
  else
    echo "PR branch and $main_branch have different hashes for submodule"\
         "$module_name:"
    echo "---- PR branch has hash $pr_hash"
    echo "---- $module_name/$main_branch has hash $main_hash"
    unsync_submodules+=( "$module_name" )
    exit_code=1
  fi
done

echo $'\n'
if [ $exit_code -eq 1 ]
then
  echo "This PR has the following submodules"\
       "that are out of sync with master or develop : ${unsync_submodules[@]}"
  echo "FAILURE : Please make sure your branch"\
       "is up to date with develop."\
       "Merge any submodule changes into the submodule's develop or"\
       "master branch before merging this PR with the main GEOSX repository."
else
  echo "SUCCESS : PR submodules are up to date!"
fi

# Renable git lfs
git lfs install &> /dev/null

exit $exit_code
