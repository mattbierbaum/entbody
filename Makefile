DIRS = examples
STUB = ${ENTBODY}

all: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

clean: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

.PHONY: all
.PHONY: clean
