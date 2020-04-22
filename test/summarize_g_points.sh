#!/bin/bash

for FILE in $@
do
    INFO=$(ncdump -h $FILE | sed -n '3,7p')
    echo $INFO $FILE
done
