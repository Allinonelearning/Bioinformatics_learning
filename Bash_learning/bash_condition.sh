#!/bin/bash

# if statement
read -p "What is your name" name
# 注意方括号内有空格空格
if [ -z "$name" ];then
	echo "Please enter your name!"
else
	echo "hello,$name"
fi

admin="Xingjinyin"
if [ "$name" == "$admin" ] ; then
	echo "you are the admin user"
else
	echo "you are not the admin user"
fi
