INSTALL_DIR = /usr/local/bin

ifdef prefix
	INSTALL_DIR = $(prefix)
endif

build:
	go build -o bin/gesyntek-run ./cmd/gesyntek-run/main.go
	go build -o bin/gesyntek-heatmap ./cmd/gesyntek-heatmap/main.go

test:
	go test -v kmer/ksplit.go kmer/ksplit_test.go

install:
	cp bin/gesyntek-run $(INSTALL_DIR)/gesyntek-run
	cp bin/gesyntek-heatmap $(INSTALL_DIR)/gesyntek-heatmap

uninstall:
	rm -f $(INSTALL_DIR)/gesyntek-run
	rm -f $(INSTALL_DIR)/gesyntek-heatmap