CPPFLAGS += -I htslib/
CFLAGS   += -g -Wall -O2  -std=c99
LDFLAGS  += $(LIBS) -lz -lm -lpthread
BUILD_DIR = build

BINARY = invar
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/invar.o \
      $(BUILD_DIR)/depth_main.o \
      $(BUILD_DIR)/thread.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/misc_p.o \
	  $(BUILD_DIR)/error.o \

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test install uninstall

$(BINARY): htslib/libhts.a $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) htslib/libhts.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/invar.o: src/invar.c src/misc.h src/error.h src/invar.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/depth_main.o: src/depth_main.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c src/invar.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc_p.o: src/misc_p.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/error.o: src/error.c src/error.h
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

test: $(BINARY)
	./scripts/test.sh

