#! /bin/bash

touch NEWS README AUTHORS Changelog

autoreconf --force --install
automake --force-missing --add-missing --copy