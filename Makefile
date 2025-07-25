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

PHONY: style
style:
	@echo "\033[0;36m[INFO] Running style formatting and lint fixes...\033[0m"
	@echo \"RUFF VERSION: `uv run ruff --version`\"
	uv run ruff check --fix-only
	uv run ruff format
	@echo "\033[0;32m[OK] Style formatting and lint completed successfully\033[0m"

PHONY: style-unsafe
style-unsafe:
	@echo "\033[0;33m[INFO] Running style formatting and lint fixes (unsafe)...\033[0m"
	@echo \"RUFF VERSION: `uv run ruff --version`\"
	uv run ruff check --fix-only --unsafe-fixes
	uv run ruff format
	@echo "\033[0;32m[OK] Style formatting and lint completed (unsafe)\033[0m"

PHONY: check
check:
	@echo "\033[0;36m[INFO] Running style checks...\033[0m"
	@echo \"RUFF VERSION: `uv run ruff --version`\"
	uv run ruff format --diff
	uv run ruff check
	@echo "\033[0;32m[OK] Style checks completed successfully\033[0m"

PHONY: test
test:
	@echo "\033[0;36m[INFO] Running tests without coverage...\033[0m"
	uv run pytest --numprocesses $(PYTEST_CORES)

PHONY: cov
cov:
	@echo "\033[0;36m[INFO] Running tests with coverage...\033[0m"
	uv run pytest --numprocesses $(PYTEST_CORES) --cov-report=term-missing --cov-config=pyproject.toml --cov=src/circuit_rl --cov=tests

PHONY: lock
lock:
	@echo "\033[0;36m[INFO] Locking dependencies...\033[0m"
	@echo \"UV VERSION: `uv --version`\"
	uv export --all-groups --no-hashes --format requirements-txt > requirements.txt
	@echo "\033[0;32m[OK] Dependencies locked successfully\033[0m"

PHONY: install
install:
	@echo "\033[0;36m[INFO] Installing project...\033[0m"
	@if ! command -v uv; then \
		echo "\033[0;31m[ERROR] uv is not installed. Please install uv first.\033[0m"; \
		echo "\033[0;33m[NOTICE] Check https://docs.astral.sh/uv/getting-started/installation/ to check how to install uv 033[0m"; \
		exit 1; \
	fi
	uv sync --all-groups
	pre-commit install
	@echo "\033[0;32m[OK] Project installed successfully\033[0m"
