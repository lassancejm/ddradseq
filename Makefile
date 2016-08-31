CC=gcc
CFLAGS=-O2
LDFLAGS=-lfasta
BARCODE_LDFLAGS=-lfasta -lz -lpthread
all: parse_pool trim_barcode trim_3prime

trim_barcode: trim_barcode.o
	$(CC) $(CFLAGS) -o $@ $^ $(BARCODE_LDFLAGS)

parse_pool: parse_pool.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

trim_3prime: trim_3prime.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f parse_pool trim_barcode trim_3prime *.o
