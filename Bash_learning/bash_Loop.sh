#!/bin/bashread -p "What is your name?" name
# for循环
user="Xingjinyin"
for user in ${user}
do
	echo ${user}
done

for num in {1..10}
do
	echo ${num}
done

# while循环
counter=1
while [ $counter -le 10 ]
do
	echo $counter
	((counter++))
done

read -p "What is your name?" name
while [ -z ${name} ]
do
	echo "Your name can not be black.Please enter your name"
	read -p "What is your name?" name
done
echo "welcome,${name}"

# until循环
count=1
until [ $count -gt 10 ]
do
	echo $count
	((count++))
done

# continue and break
for i in 1 2 3 4 5
do
	if [ $i -eq 2 ]
	then
		echo "skipping number 2"
		continue
	fi
	echo "i is equal to $i"
done

num=1
while [ $num -lt 10 ]
do
	if [ $num -eq 5 ]
	then
		break
	fi
	((num++))
done
echo "LOOP completed"













