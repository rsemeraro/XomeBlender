PICARD_CHECK = $(shell find / -type f -name 'picard.jar' -print -quit > deps.txt 2>/dev/null)
GATK_CHECK = $(shell find /  -type f -name 'GenomeAnalysisTK.jar' -print -quit >> deps.txt 2>/dev/null)
VCFTOOLS_CHECK = $(shell find / -type f -name 'vcftools' -print -quit >> deps.txt 2>/dev/null)
GCOMP_VERSION = $(shell g++ -dumpversion | cut -f1-2 -d.)
GC_4_7 = $(shell echo $(GCOMP_VERSION)\>4.6 | bc )
ARCH= -static
ARCHECK = $(shell uname -s)
STD_VAR = c++0x
ifeq ($(GC_4_7),1)
STD_VAR = c++11
endif
ifeq ($(ARCHECK), Darwin)
ARCH =
endif
all: src/inx_counter.cc src/xome_counter.cc
	g++ $(ARCH) -march=native -mtune=native -O2 -std=$(STD_VAR) -fomit-frame-pointer -Wextra src/inx_counter.cc -o Scripts/inx_counter
	g++ $(ARCH) -march=native -mtune=native -O2 -std=$(STD_VAR) -fomit-frame-pointer -Wextra src/xome_counter.cc -o Scripts/xome_counter
	g++ $(ARCH) -march=native -mtune=native -O2 -std=$(STD_VAR) -fomit-frame-pointer -Wextra src/count_vars.cc -o Scripts/count_vars
	g++ $(ARCH) -march=native -mtune=native -O2 -std=$(STD_VAR) -fomit-frame-pointer -Wextra src/randlines.cc -o Scripts/randlines
	$(PICARD_CHECK)
	$(GATK_CHECK)
	$(VCFTOOLS_CHECK)	
clean: 
	$(RM) alt_counter
