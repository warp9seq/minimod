VERSION = $(shell git describe --tags --dirty 2>/dev/null || echo "0.0.0-unknown")
VERSION = v0.2.0

CPPFLAGS += -I htslib/
CFLAGS   += -g -Wall -O2  -std=c99 -DMINIMOD_VERSION=\"$(VERSION)\"
LDFLAGS  += $(LIBS) -lz -lm -lpthread
BUILD_DIR = build

BINARY = minimod
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/minimod.o \
      $(BUILD_DIR)/view_main.o \
	  $(BUILD_DIR)/mod_freq_main.o \
      $(BUILD_DIR)/thread.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/misc_p.o \
	  $(BUILD_DIR)/error.o \
	  $(BUILD_DIR)/mod.o \
	  $(BUILD_DIR)/ref.o

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test install uninstall

$(BINARY): htslib/libhts.a $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) htslib/libhts.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c src/minimod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/minimod.o: src/minimod.c src/misc.h src/error.h src/minimod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/view_main.o: src/view_main.c src/error.h src/minimod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/mod_freq_main.o: src/mod_freq_main.c src/error.h src/minimod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c src/minimod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc_p.o: src/misc_p.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/error.o: src/error.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/mod.o: src/mod.c src/mod.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/ref.o: src/ref.c src/kseq.h src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

htslib/libhts.a:
	@if test -e $(BUILD_DIR)/lib/libhts.a; then \
		echo "htslib found at htslib/libhts.a"; \
	else \
		echo "htslib not found at htslib/libhts.a"; \
		echo "Please run 'scripts/install-hts.sh' first"; \
		exit 1; \
	fi

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/*

release: distclean
# make the release
	mkdir -p minimod-$(VERSION)/build
	cp -r README.md LICENSE Makefile scripts src docs minimod-$(VERSION)
	tar -zcf minimod-$(VERSION)-release.tar.gz minimod-$(VERSION)
	rm -rf minimod-$(VERSION)
# make the binaries
	make -j8
	mkdir -p minimod-$(VERSION)
	mv minimod minimod-$(VERSION)/build
	cp -r README.md LICENSE minimod-$(VERSION)/
	tar -zcf minimod-$(VERSION)-x86_64-linux-binaries.tar.gz minimod-$(VERSION)
	rm -rf minimod-$(VERSION)
	tar xf minimod-$(VERSION)-x86_64-linux-binaries.tar.gz
	mv minimod-$(VERSION)/minimod minimod
	rm -rf minimod-$(VERSION)
	test/test.sh

test: $(BINARY)
	./test/test.sh

memtest: $(BINARY)
	./test/test.sh mem

