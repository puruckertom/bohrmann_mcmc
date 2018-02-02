#!/bin/csh
free -m | grep "^Mem" | awk '{print $2}'
