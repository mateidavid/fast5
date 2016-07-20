.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all help clean check_virtualenv develop develop-user develop-uninstall develop-uninstall-user

PYTHON = $(shell which python)

all: help

print-%:
	@echo '$*=$($*)'

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

clean: ## Remove build products
	${PYTHON} setup.py clean
	rm -rf fast5.egg-info build dist
	find fast5/ \( -name '*.pyc' -o -name '*.so' \) -delete

check_virtualenv:
	@[ "$$VIRTUAL_ENV" ] || { echo "not in a virtualenv" >&2; exit 1; }

develop: check_virtualenv clean ## Install in develop mode to current virtualenv
	${PYTHON} setup.py develop

develop-user: clean ## Install in develop mode to current user
	${PYTHON} setup.py develop --user

develop-uninstall: check_virtualenv clean ## Uninstall from current virtualenv
	${PYTHON} setup.py develop --uninstall

develop-uninstall-user: clean ## Uninstall from current user
	${PYTHON} setup.py develop --uninstall --user
