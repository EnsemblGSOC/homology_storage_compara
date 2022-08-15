#!/bin/sh
for FILE in $1/*; do
    filename=$(basename $FILE)
    bin/gtutils -s $2/$filename $FILE
done