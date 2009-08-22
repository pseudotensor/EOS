#!/bin/sh

echo "First"
echo $1
echo $2

echo "$1"
echo "$2"

#sh testarg2.sh boob "$*"
sh testarg2.sh boob "$@"
