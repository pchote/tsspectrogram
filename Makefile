# Copyright 2010, 2011, 2012 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

USE_READLINE := TRUE

CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_POSIX_C_SOURCE=200112L -D_BSD_SOURCE -I/usr/local/include/
LFLAGS = -lcpgplot -lpgplot -lm -L/usr/local/lib/

ifeq ($(shell uname),Darwin)
	LINKER = gcc
	CFLAGS += -D_DARWIN_C_SOURCE
endif

SRC = tsspectrogram.c
OBJ = $(SRC:.c=.o)

tsspectrogram: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) tsspectrogram

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
