#-----------------------------#
PROJECT(GASCOIGNESTD)
#-----------------------------#

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#
# .gascoignerc contains the user dependent settings
#

INCLUDE($ENV{HOME}/gascoignerc_example)

#
# CMakeGlobal.txt is intended for use in derived projects
#

INCLUDE(CMakeGlobal.txt)

#
# the following is local
#

SUBDIRS(src test Examples doc)
