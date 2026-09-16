#CFLAGS=     -g -Wall -fsanitize=address -fno-omit-frame-pointer -fno-strict-aliasing -Wno-unused-function -Wno-deprecated-declarations -Wno-array-bounds
CFLAGS=		-O3 -Wall -fno-strict-aliasing -Wno-unused-function -Wno-deprecated-declarations -Wno-array-bounds
CPPFLAGS=
INCLUDES=
LDFLAGS=
OBJS=
PROG=       hapscotch hapcount hapcure hictools seqtools
PROG_EXTRA=
LIBS=       -lm -lz -lpthread
DESTDIR=    ~/bin

HiGHS_ROOT = HiGHS
HiGHS_DIR = $(HiGHS_ROOT)/build
LIBHiGHS_VERSION_MAJOR = 1
LIBHiGHS_VERSION_MINOR = 15

LIB_DIR = lib

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SHARED_EXT = dylib
	LIBHiGHS_EXT_MAJOR = $(LIBHiGHS_VERSION_MAJOR).$(SHARED_EXT)
	LIBHiGHS_EXT_MINOR = $(LIBHiGHS_VERSION_MAJOR).$(LIBHiGHS_VERSION_MINOR).$(SHARED_EXT)
	R_PATH = -L$(LIB_DIR) -Wl,-rpath,@executable_path/$(LIB_DIR)
	MACOSX_DEPLOYMENT_TARGET ?= $(shell sw_vers -productVersion | cut -d. -f1-2)
	CFLAGS += -mmacosx-version-min=$(MACOSX_DEPLOYMENT_TARGET)
	CMAKE_FLAGS += -DCMAKE_OSX_DEPLOYMENT_TARGET=$(MACOSX_DEPLOYMENT_TARGET)
else
	SHARED_EXT = so
	LIBHiGHS_EXT_MAJOR = $(SHARED_EXT).$(LIBHiGHS_VERSION_MAJOR)
	LIBHiGHS_EXT_MINOR = $(SHARED_EXT).$(LIBHiGHS_VERSION_MAJOR).$(LIBHiGHS_VERSION_MINOR)
	R_PATH = -L$(LIB_DIR) -Wl,-rpath,'$$ORIGIN/$(LIB_DIR)'
endif

HiGHS_INCLUDES = -I$(HiGHS_ROOT) -I$(HiGHS_ROOT)/highs -I$(HiGHS_DIR)
HiGHS_LIB = $(HiGHS_DIR)/lib/libhighs.$(SHARED_EXT)
HiGHS_LIBS = -lhighs

HiGHS_OBJS = $(LIB_DIR)/libhighs.$(SHARED_EXT) \
	$(LIB_DIR)/libhighs.$(LIBHiGHS_EXT_MINOR) \
	$(LIB_DIR)/libhighs.$(LIBHiGHS_EXT_MAJOR)

.PHONY: all extra clean depend test highs

all: $(HiGHS_OBJS) $(PROG)

extra: all $(PROG_EXTRA)

debug: $(PROG)
debug: CFLAGS += -DDEBUG

%.o: %.c | $(HiGHS_LIB)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(HiGHS_INCLUDES) $< -o $@

highs: $(HiGHS_LIB)

$(HiGHS_LIB):
	+cd $(HiGHS_ROOT) && cmake -S . -B build $(CMAKE_FLAGS) && \
	cmake --build build --target highs

$(LIB_DIR):
	@mkdir -p $@

# Copy HiGHS library to lib directory
$(LIB_DIR)/libhighs.$(SHARED_EXT): $(HiGHS_LIB) | $(LIB_DIR)
	@cp -f $< $@
ifeq ($(UNAME_S),Darwin)
	@install_name_tool -id @rpath/libhighs.$(SHARED_EXT) $@
endif

$(LIB_DIR)/libhighs.$(LIBHiGHS_EXT_MAJOR): $(LIB_DIR)/libhighs.$(SHARED_EXT)
	@ln -sf libhighs.$(SHARED_EXT) $@

$(LIB_DIR)/libhighs.$(LIBHiGHS_EXT_MINOR): $(LIB_DIR)/libhighs.$(LIBHiGHS_EXT_MAJOR)
	@ln -sf libhighs.$(LIBHiGHS_EXT_MAJOR) $@

hapscotch: hapscotch.o alnio.o busco.o overlap.o ploidy.o hap.o hic.o sdict.o paf.o range.o cov.o bamlite.o ONElib.o misc.o kthread.o kalloc.o kopen.o | $(HiGHS_OBJS)
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -L. $(R_PATH) $(HiGHS_LIBS) $(LIBS)

hapcount: hapcount.o alnio.o busco.o overlap.o ploidy.o sdict.o paf.o range.o misc.o kthread.o kalloc.o kopen.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ -L. $(LIBS)

hapcure: hapcure.o ec.o hic.o sdict.o cov.o bamlite.o ONElib.o misc.o kalloc.o kopen.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ -L. $(LIBS)

hictools: hictools.o hic.o sdict.o cov.o bamlite.o ONElib.o misc.o kalloc.o kopen.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ -L. $(LIBS)

seqtools: seqtools.o sdict.o bgzf.o misc.o kalloc.o kopen.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ -L. $(LIBS)

clean:
	rm -fr *.o a.out $(PROG) $(OBJS) $(PROG_EXTRA)

install: | $(HiGHS_LIB)
	cp $(PROG) $(DESTDIR)
	cp -Paf $(LIB_DIR) $(DESTDIR)

depend: | $(HiGHS_LIB)
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

sdict.o: sdict.h agp-spec.h misc.h khash.h ksort.h kseq.h kvec.h
paf.o: paf.h misc.h kseq.h
cov.o: cov.h
bamlite.o: bamlite.h
bgzf.o: bgzf.h
range.o: range.h misc.h kavl.h
misc.o: misc.h kseq.h
kthread.o: kthread.h
kalloc.o: kalloc.h
ONElib.o: ONElib.h
hapscotch.o: paf.h hic.h sdict.h agp-spec.h misc.h ketopt.h kvec.h kseq.h khash.h kthread.h busco.h overlap.h ploidy.h hap.h alnio.h version.h
hapcount.o: paf.h sdict.h misc.h ketopt.h kvec.h busco.h overlap.h ploidy.h alnio.h version.h
alnio.o: alnio.h misc.h sdict.h overlap.h paf.h kvec.h
busco.o: busco.h khash.h sdict.h agp-spec.h misc.h kseq.h kvec.h
overlap.o: overlap.h misc.h kvec.h kthread.h sdict.h
ploidy.o: ploidy.h overlap.h misc.h kvec.h sdict.h range.h
hap.o: hap.h busco.h overlap.h ploidy.h misc.h sdict.h agp-spec.h range.h hic.h kvec.h kthread.h
hic.o: hic.h sdict.h cov.h misc.h ketopt.h kvec.h kseq.h khash.h ONElib.h
ec.o: ec.h hic.h sdict.h misc.h kvec.h
hictools.o: hic.h sdict.h misc.h bamlite.h khash.h kstring.h kvec.h ketopt.h version.h
hapcure.o: ec.h hic.h sdict.h misc.h ketopt.h version.h
seqtools.o: agp-spec.h sdict.h bgzf.h misc.h kvec.h kstring.h ketopt.h version.h
