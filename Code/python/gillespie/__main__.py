# __main__.py
import os, sys

from importlib import import_module


def main():
    # dynamic dispatch, e.g. based on sys.argv
    name = sys.argv[1]
    mod = import_module(f"mypkg.{name}")
    getattr(mod, "main")()


if __name__ == "__main__":
    main()
