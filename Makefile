INSTALL_DIR = /usr/local/bin

ifdef prefix
	INSTALL_DIR = $(prefix)
endif

build:
	go build -o bin/gesyntek-run ./cmd/gesyntek-run/main.go

test:
	go test -v kmer/ksplit.go kmer/ksplit_test.go