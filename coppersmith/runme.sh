#!/bin/bash

DIRECTORY=$(cd `dirname $0` && pwd);
gcc -o qrfactor Householder_QR.c -std=c99 -I$DIRECTORY/.. -L/$DIRECTORY/.. -lgmp -lmpfr -lflint -Wl,-rpath,$DIRECTORY/..;
gcc flint_copper.c -o copper_prog -std=c99 -I$DIRECTORY/.. -L/$DIRECTORY/.. -lmpfr -lgmp -lflint -Wl,-rpath,$DIRECTORY/..;
