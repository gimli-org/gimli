#! /bin/bash

g++ -shared -fPIC `wx-config --cflags --libs` `pkg-config gtk+-2.0 --cflags --libs` -lGL getHandleHack.cpp -o libgetHandleHack.so

