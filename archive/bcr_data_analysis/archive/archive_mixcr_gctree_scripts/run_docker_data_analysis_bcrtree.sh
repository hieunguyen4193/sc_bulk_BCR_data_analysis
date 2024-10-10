#!/bin/bash
port=$1;
name=$2;
password=$3
sudo docker run -it --cpus 20 --memory 64GB \
-v /home/hieunguyen:/home/hieunguyen -v /media:/media -v /mnt:/mnt -p ${port}:8787 -e PASSWORD=${password} --name ${name} \
tronghieunguyen/bcrtree:latest
