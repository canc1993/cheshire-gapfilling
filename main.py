from predict import *
from validate import *
import logging

logging.disable()


def main():
    # merge reaction pool with GEMs
    create_pool()
    # predict scores for reactions in reaction pool
    predict()
    # predict phenotypes
    validate()


if __name__ == "__main__":
    main()