#!/bin/bash

# This script is used to get the number of line changes with time
# Usage:
#   bash repositoryStats.sh
#

months=`seq 1 12`
years=`seq 2016 2024`

echo "Month-Year NumDevelopersMonth numDevelopers3MonthWindow NumDevelopersTotal NumCommitsMonth NumCommitsTotal NumLinesInMonth"
for year in $years
do
    for month in $months
    do
      numLines=$(git log --since="$year-$month-01" --until="$year-$month-31" --format= --numstat | awk '{s+=$1; s+=$2} END {print s}')
      #numDevelopersMonth=$(git log --since="$year-$month-01" --until="$year-$month-31" --all --pretty="%an" | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
      numDevelopersMonth=$(git log --since="$year-$month-01" --until="$year-$month-31" | grep -E 'Author:|Co-authored-by:' | sed 's/^.*: //' | sed 's/<.*//' | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
      
      minus1month=$(date -j -v-1m -f "%Y-%m-%d" "$year-$month-15" "+%Y-%m-%d")
      plus1month=$(date -j -v+1m -f "%Y-%m-%d" "$year-$month-15" "+%Y-%m-%d")

#      echo "Begin: $minus1month End: $plus1month"

      # numDevelopers3Month=$(git log --since="$minus1month"  --until="$plus1month" --all --pretty="%an" | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
      # numDevelopersTotal=$(git log --since="2010-01-01"      --until="$year-$month-31" --all --pretty="%an" | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
      numDevelopers3Month=$(git log --since="$minus1month"  --until="$plus1month" | grep -E 'Author:|Co-authored-by:' | sed 's/^.*: //' | sed 's/<.*//' | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
      numDevelopersTotal=$(git log --since="2010-01-01"      --until="$year-$month-31" | grep -E 'Author:|Co-authored-by:' | sed 's/^.*: //' | sed 's/<.*//' | sort | uniq | tr -d " ' " | xargs python3 repositoryStats_condenseDev.py)
 
      numCommitsMonth=$(git rev-list --count HEAD --since="$year-$month-01"  --before="$year-$month-31")
      numCommitsTotal=$(git rev-list --count HEAD --since="2010-01-01"  --before="$year-$month-31")
      echo "$month-$year $numDevelopersMonth $numDevelopers3Month $numDevelopersTotal $numCommitsMonth $numCommitsTotal $numLines"
    done
done


