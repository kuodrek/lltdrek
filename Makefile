.PHONY: install test lint

install:
	poetry install

test:
	poetry run pytest

lint:
	poetry run pre-commit run --all-files
