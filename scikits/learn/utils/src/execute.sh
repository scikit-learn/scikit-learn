#! /bin/sh
WEKAJAR=$PWD/weka.jar
CLASSPATH=$WEKAJAR:$PWD java testarff $1
