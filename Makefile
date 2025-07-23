PHONY: build-docs
build-docs:
	#cp AUTHORS.md docs/source/project/
	#cp CHANGELOG.md docs/source/project/
	#cp CONTRIBUTING.md docs/source/project/
	cp README.md docs/source/project/
	uv run jupyter-book config sphinx docs/
	uv run jupyter-book build docs/
	uv pip freeze > docs/requirements.txt

PHONY: serve
serve:
	uv run python3 -m http.server --directory docs/_build/html