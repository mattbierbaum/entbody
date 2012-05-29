DIRS = brazilnuts polycrystals polycrystals_morse rayleightaylor sand stir moshpits springs source_sink
STUB = ${ENTBODY}/examples

all: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

clean: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

.PHONY: all
.PHONY: clean
