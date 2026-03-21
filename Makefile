.PHONY: install test lint

install:
	pip install -e ".[dev]"

test:
	pytest

lint:
	pre-commit run --all-files
