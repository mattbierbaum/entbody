DIRS = `find . -maxdepth 1 -mindepth 1 -type d` 
STUB = ${ENTBODY}/examples

all: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

clean: 
	for dir in $(DIRS);  do cd "$$dir" && $(MAKE) $@ && cd ..; done

.PHONY: all
.PHONY: clean
