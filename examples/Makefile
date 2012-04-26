DIRS = brazilnuts polycrystals rayleightaylor sand stir
STUB = ${ENTBODY}/examples

all: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

clean: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

.PHONY: all
.PHONY: clean
