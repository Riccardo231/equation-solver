#!/bin/bash

while IFS='|' read -r equation expected; do
  echo "Testing: $equation"
  output=$(./main "$equation")
  if [ -n "$output" ]; then
    echo "PASS: got '$output'"
  else
    echo "FAIL: no output"
  fi
done < testcases.txt
