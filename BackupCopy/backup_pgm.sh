#!/bin/sh

folder='/media/PGM/results_S5/'

if [ $1 = 'cp' ]
then
	cp "$2" "$folder$3"
fi

if [ $1 = 'mkdir' ]
then
	mkdir "$folder$2"
fi

exit 0
