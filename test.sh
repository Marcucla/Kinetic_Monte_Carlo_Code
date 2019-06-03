#!/bin/bash
sudo rm -rf msr
sudo mkdir msr
cd ./msr
sudo cp ../test3new1.config ./
sudo cmake ../src
sudo make -j4
sudo ./DislocationKMC --config test3new1.config
