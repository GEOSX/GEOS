#!/bin/bash

DIFF=$(git diff --name-only origin/develop)

if [[ $DIFF == *".integrated_tests.yaml"* ]]; then
  echo "Changes to the integrated test baseline detected"
  if [[ $DIFF != *"BASELINE_NOTES.md"* ]]; then
    echo "You need to add a note to the BASELINE_NOTES.md file to justify changes to the test baselines"
    exit 1
  fi
fi

