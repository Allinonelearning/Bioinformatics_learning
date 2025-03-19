#!/bin/bash
# This script displays the date and who's logged on
date
who
# copy the /usr/bin directory listing to a log file
today=$(date +%y%m%d)
ls /usr/bin -al > log.$today
