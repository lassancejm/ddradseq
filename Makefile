# Makefile for the ddradseq program
# Lummei Analytics LLC November 2016
CC = gcc
CFLAGS = -O2 -Wall -D_FILE_OFFSET_BITS=64 -std=gnu99
DEBUG_CFLAGS = -ggdb -Wall -D_FILE_OFFSET_BITS=64 -std=gnu99
LDFLAGS = -lz
TARGET = ddradseq
SRCS = $(wildcard *.c)
DEPS = $(wildcard *.h)
OBJS = $(SRCS:.c=.o)
INSTALL_DIR = $(HOME)/bin

all: $(TARGET)

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)
	@echo "Done"

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean install

install:
	cp $(TARGET) $(INSTALL_DIR)

clean:
	rm -f $(TARGET) $(OBJS)
