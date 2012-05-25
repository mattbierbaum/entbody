DIRS = brazilnuts polycrystals polycrystals_morse rayleightaylor sand stir moshpits springs
STUB = ${ENTBODY}/examples

all: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

clean: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

.PHONY: all
.PHONY: clean
