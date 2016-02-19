#!/bin/bash

export CFLAGS="-I$PREFIX/include"
$PYTHON setup.py install
